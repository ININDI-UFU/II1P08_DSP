// Soma de 3 senoides (250 Hz + 800 Hz + 1.5 kHz) amostradas a fs=8 kHz,
// FFT do sinal resultante (esp-dsp) e visualização via wserial (LasecPlot/Teleplot).
#include <Arduino.h>
#include <math.h>

#include "services/wserial.h"
#include "dsps_fft2r.h"
#include "dsps_tone_gen.h"
#include "dsps_add.h"
#include "dsps_view.h"
#include "dsps_biquad.h"
#include "dsps_biquad_gen.h"
#include "dsps_wind_hann.h"

// ---------- Parâmetros do sinal ----------
constexpr int       N  = 1024;           // número de amostras (potência de 2 -> exigido pela FFT real)
constexpr float     fs = 8000.0f;        // frequência de amostragem [Hz]
constexpr float     f1 = 250.0f;         // tom 1 [Hz]
constexpr float     f2 = 800.0f;         // tom 2 [Hz]
constexpr float     f3 = 1500.0f;        // tom 3 [Hz]

// dsps_tone_gen_f32 usa frequência normalizada pela Nyquist (fs/2): freq_norm = f / (fs/2)
constexpr float toNyquist(float f) { return 2.0f * f / fs; }

static float tone1[N];
static float tone2[N];
static float tone3[N];
static float real_signal[N]; // soma das 3 senoides (sinal real)

// Buffer de trabalho da FFT: pares [re0, im0, re1, im1, ...] -> tamanho 2*N
__attribute__((aligned(16))) static float fft_buf[N * 2];

void setup() {
    wserial.begin(115200);
    delay(1000);

    // Inicializa as tabelas de coeficientes da FFT (uma vez só)
    esp_err_t ret = dsps_fft2r_init_fc32(NULL, N);
    if (ret != ESP_OK) {
        wserial.println(String("Falha ao inicializar FFT: ") + ret);
        return;
    }

    // ---------- Gera as 3 senoides e soma em um único sinal real ----------
    dsps_tone_gen_f32(tone1, N, 1.0f, toNyquist(f1), 0.0f);
    dsps_tone_gen_f32(tone2, N, 1.0f, toNyquist(f2), 0.0f);
    dsps_tone_gen_f32(tone3, N, 1.0f, toNyquist(f3), 0.0f);

    dsps_add_f32(tone1, tone2, real_signal, N, 1, 1, 1); // real_signal = tone1 + tone2
    dsps_add_f32(real_signal, tone3, real_signal, N, 1, 1, 1); // real_signal += tone3

    // ---------- Filtro passa-faixa: isola apenas o tom do meio (f2 = 800 Hz) ----------
    // dsps_biquad_gen_* espera a frequência normalizada por fs (0..0.5), diferente do
    // toNyquist() usado no tone_gen (que normaliza por fs/2) — aqui é direto f2/fs.
    // "bpf0db" mantém ganho 0 dB na frequência central (não atenua o tom desejado).
    // Cascata de 2 biquads (cada um com seu próprio estado w[]) para rejeitar melhor
    // os tons vizinhos (250 Hz e 1.5 kHz), que ficam fora da banda passante mas
    // seriam pouco atenuados por um único estágio de 2ª ordem.
    constexpr float bpfQ = 5.0f;
    static float bpfCoeffs[5];
    dsps_biquad_gen_bpf0db_f32(bpfCoeffs, f2 / fs, bpfQ);

    static float filtered_signal[N]; // sinal filtrado: deve restar só o tom de 800 Hz
    static float bpf_stage1[N];
    float bpf_w1[2] = {0.0f, 0.0f};
    float bpf_w2[2] = {0.0f, 0.0f};
    dsps_biquad_f32(real_signal, bpf_stage1, N, bpfCoeffs, bpf_w1);
    dsps_biquad_f32(bpf_stage1, filtered_signal, N, bpfCoeffs, bpf_w2);

    // ---------- Janela (reduz o vazamento espectral antes da FFT) ----------
    // Sem janela, a FFT assume implicitamente uma janela retangular: tons que não
    // caem exatamente em um bin (como f2=800 Hz, já que 800/df=102.4 não é inteiro)
    // "vazam" para os bins vizinhos. A janela de Hann suaviza as bordas do bloco de
    // N amostras e reduz bastante esse vazamento (ao custo de um lóbulo principal
    // um pouco mais largo).
    static float hann[N];
    static float windowed_signal[N];
    dsps_wind_hann_f32(hann, N);
    for (int i = 0; i < N; i++) {
        windowed_signal[i] = real_signal[i] * hann[i];
    }

    // esp-dsp não tem uma função para empacotar um vetor real no array complexo
    // intercalado da FFT — esse passo é sempre manual: parte imaginária zerada.
    for (int i = 0; i < N; i++) {
        fft_buf[2 * i + 0] = windowed_signal[i]; // real (já com janela aplicada)
        fft_buf[2 * i + 1] = 0.0f;               // imaginário
    }

    // ---------- FFT ----------
    dsps_fft2r_fc32(fft_buf, N);              // FFT radix-2 in-place
    dsps_bit_rev_fc32(fft_buf, N);             // reordena bit-reversal
    dsps_cplx2reC_fc32(fft_buf, N);            // converte para espectro real (0..N/2)

    // ---------- Magnitude em dB ----------
    static float magnitude[N / 2];
    for (int i = 0; i < N / 2; i++) {
        const float re = fft_buf[2 * i + 0];
        const float im = fft_buf[2 * i + 1];
        magnitude[i] = 10 * log10f((re * re + im * im) / (N * N) + 1e-12f);
    }

    // ---------- Visualização via wserial (LasecPlot/Teleplot) ----------
    wserial.println("Espectro (FFT) do sinal 250 Hz + 800 Hz + 1.5 kHz @ fs=8 kHz");

    // Sinal no tempo: dt_ms=1 aproxima o eixo X por índice de amostra
    // (o período real de amostragem, 1/fs = 0.125 ms, não é representável em ms inteiro).
    wserial.plot("signal", 1, real_signal, N, "V");
    wserial.plot("filtered", 1, filtered_signal, N, "V");

    // Espectro de magnitude: dt_ms=1 aproxima o eixo X por índice de bin da FFT.
    wserial.plot("magnitude", 1, magnitude, N / 2, "dB");

    // Gráfico ASCII do espectro.
    dsps_view(magnitude, N / 2, 128, 20, -80, 0, '*');
}

void loop() {
    wserial.update();
}

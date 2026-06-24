// Soma de 3 senoides (250 Hz + 800 Hz + 1.5 kHz) amostradas a fs=8 kHz,
// filtro passa-faixa, FFT do sinal resultante (esp-dsp) e visualização via
// wserial (LasecPlot/Teleplot). Cada etapa do processamento DSP é uma função
// própria, para ficar fácil de seguir e explicar isoladamente.
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
constexpr int   N  = 1024;    // número de amostras (potência de 2 -> exigido pela FFT real)
constexpr float fs = 8000.0f; // frequência de amostragem [Hz]
constexpr float f1 = 250.0f;  // tom 1 [Hz]
constexpr float f2 = 800.0f;  // tom 2 [Hz] -> é o tom que o filtro passa-faixa isola
constexpr float f3 = 1500.0f; // tom 3 [Hz]

// dsps_tone_gen_f32 usa frequência normalizada pela Nyquist (fs/2): freq_norm = f / (fs/2)
constexpr float toNyquist(float f) { return 2.0f * f / fs; }

// Buffers compartilhados entre as etapas (preenchidos em setup() e enviados ao final)
static float real_signal[N];          // soma dos 3 tons, no tempo
static float filtered_signal[N];      // saída do filtro passa-faixa, no tempo (deve restar só o tom de 800 Hz)
static float magnitude[N / 2];        // espectro (frequência) do sinal somado, em dB
static float magnitude_filtered[N / 2]; // espectro (frequência) do sinal filtrado, em dB

// ---------- Etapa 1: gerar o sinal de teste (3 tons somados) ----------
void gerarSinal(float *saida) {
    static float tone1[N];
    static float tone2[N];
    static float tone3[N];

    dsps_tone_gen_f32(tone1, N, 1.0f, toNyquist(f1), 0.0f);
    dsps_tone_gen_f32(tone2, N, 1.0f, toNyquist(f2), 0.0f);
    dsps_tone_gen_f32(tone3, N, 1.0f, toNyquist(f3), 0.0f);

    dsps_add_f32(tone1, tone2, saida, N, 1, 1, 1); // saida = tone1 + tone2
    dsps_add_f32(saida, tone3, saida, N, 1, 1, 1); // saida += tone3
}

// ---------- Etapa 2: filtro passa-faixa (isola só o tom central, f2) ----------
// Cascata de 2 biquads IIR idênticos (cada um com seu próprio estado w[]): a
// resposta em frequência fica ao quadrado, dobrando a atenuação (em dB) fora
// da banda passante sem perder ganho no centro ("bpf0db" = 0 dB em f2).
void filtrarPassaFaixa(const float *entrada, float *saida) {
    constexpr float bpfQ = 5.0f;
    static float bpfCoeffs[5];
    // dsps_biquad_gen_* normaliza a frequência por fs direto (0..0.5) -- diferente
    // do toNyquist() usado em gerarSinal(), que normaliza por fs/2.
    dsps_biquad_gen_bpf0db_f32(bpfCoeffs, f2 / fs, bpfQ);

    static float estagio1[N];
    float w1[2] = {0.0f, 0.0f};
    float w2[2] = {0.0f, 0.0f};
    dsps_biquad_f32(entrada, estagio1, N, bpfCoeffs, w1); // 1º estágio
    dsps_biquad_f32(estagio1, saida, N, bpfCoeffs, w2);    // 2º estágio (em série)
}

// ---------- Etapa 3: espectro de magnitude em dB (janela + FFT) ----------
// Sem janela, a FFT assume implicitamente uma janela retangular: tons que não
// caem exatamente em um bin (como f2=800 Hz, já que 800/df=102.4 não é inteiro)
// "vazam" para os bins vizinhos. A janela de Hann reduz bastante esse vazamento.
void calcularEspectroDb(const float *sinal, float *magnitudeDb) {
    // Buffer de trabalho da FFT: pares [re0, im0, re1, im1, ...] -> tamanho 2*N
    static float fft_buf[N * 2];
    static float hann[N];

    // dsps_mul_f32/dsps_memset com stride != 1 (escrever direto nas posições
    // pares do buffer intercalado) chegaram a corromper memória vizinha nessa
    // versão do esp-dsp -- a doc já avisa "step deveria ser 1 por padrão", ou
    // seja, stride != 1 não é garantido na variante otimizada em assembly.
    // Por segurança, janela e empacotamento ficam em loop manual (stride 1 só).
    dsps_wind_hann_f32(hann, N);
    for (int i = 0; i < N; i++) {
        fft_buf[2 * i + 0] = sinal[i] * hann[i]; // real, já com janela aplicada
        fft_buf[2 * i + 1] = 0.0f;               // imaginário
    }

    dsps_fft2r_fc32(fft_buf, N);   // FFT radix-2 in-place
    dsps_bit_rev_fc32(fft_buf, N); // reordena bit-reversal
    dsps_cplx2reC_fc32(fft_buf, N); // converte para espectro real (0..N/2)

    // Mesmo motivo acima: magnitude (re²+im²) também em loop manual, sem stride.
    for (int i = 0; i < N / 2; i++) {
        const float re = fft_buf[2 * i + 0];
        const float im = fft_buf[2 * i + 1];
        magnitudeDb[i] = 10 * log10f((re * re + im * im) / (N * N) + 1e-12f);
    }
}

// ---------- Etapa 4: enviar tudo para visualização ----------
// 4 séries: sinal somado e sinal filtrado no tempo, e os respectivos espectros
// (magnitude em dB) na frequência -- para comparar diretamente antes/depois do filtro.
void plotarResultados(const float *sinal, const float *filtrado,
                       const float *magnitudeDb, const float *magnitudeFiltradoDb) {
    wserial.println("Espectro (FFT) do sinal 250 Hz + 800 Hz + 1.5 kHz @ fs=8 kHz");

    // dt_ms=1 aproxima o eixo X por índice de amostra/bin (o período real de
    // amostragem, 1/fs = 0.125 ms, não é representável em ms inteiro).
    wserial.plot("signal", 1, sinal, N, "V");             // tempo: sinal somado
    wserial.plot("filtered", 1, filtrado, N, "V");         // tempo: sinal filtrado
    wserial.plot("magnitude", 1, magnitudeDb, N / 2, "dB");                 // frequência: sinal somado
    wserial.plot("magnitude_filtered", 1, magnitudeFiltradoDb, N / 2, "dB"); // frequência: sinal filtrado

    dsps_view(magnitudeDb, N / 2, 128, 20, -80, 0, '*'); // gráfico ASCII do espectro (sinal somado)
}

void setup() {
    wserial.begin(115200);
    delay(1000);

    esp_err_t ret = dsps_fft2r_init_fc32(NULL, N); // inicializa as tabelas de coeficientes da FFT (uma vez só)
    if (ret != ESP_OK) {
        wserial.println(String("Falha ao inicializar FFT: ") + ret);
        return;
    }

    gerarSinal(real_signal);
    filtrarPassaFaixa(real_signal, filtered_signal);
    calcularEspectroDb(real_signal, magnitude);
    calcularEspectroDb(filtered_signal, magnitude_filtered);
    plotarResultados(real_signal, filtered_signal, magnitude, magnitude_filtered);
}

void loop() {
    wserial.update();
}

// Soma de 3 senoides (250 Hz + 800 Hz + 1.5 kHz) amostradas a fs=8 kHz,
// filtro passa-faixa, FFT do sinal resultante (esp-dsp) e visualizacao via
// wserial (LasecPlot/Teleplot). Cada etapa do processamento DSP e uma funcao
// propria, para ficar facil de seguir e explicar isoladamente.
#include <Arduino.h>
#include <math.h>

#include "services/wserial.h"
#include "dsps_fft2r.h"
#include "dsps_tone_gen.h"
#include "dsps_add.h"
#include "dsps_view.h"
#include "dsps_biquad.h"
#include "dsps_wind_hann.h"

// ---------- Parametros do sinal ----------
constexpr int   N  = 1024;    // numero de amostras (potencia de 2 -> exigido pela FFT real)
constexpr float fs = 8000.0f; // frequencia de amostragem 8K[Hz]
constexpr float f1 = 100.0f;  // tom 1 [Hz] (100 Hz)
constexpr float f2 = 800.0f;  // tom 2 [Hz] -> e o tom que o filtro passa-faixa isola
constexpr float f3 = 2200.0f; // tom 3 [Hz] (mais distante da banda de 800 Hz)

// dsps_tone_gen_f32 usa freq_norm = f/fs (confirmado empiricamente: o periodo
// amostrado bate com f/fs, nao com f/(fs/2) -- a doc da esp-dsp chama o
// parametro de "Nyquist frequency -1..1", o que sugeria f/(fs/2), mas isso
// estava dobrando a frequencia de todo tom gerado, ex.: "800 Hz" saia a 1600 Hz).
constexpr float toFreqNorm(float f) { return f / fs; }

// Buffers compartilhados entre as etapas (preenchidos em setup() e enviados ao final)
static float real_signal[N];            // soma dos 3 tons, no tempo
static float filtered_signal[N];        // saida do filtro passa-faixa, no tempo
static float magnitude[N / 2];          // espectro do sinal somado, em dB
static float magnitude_filtered[N / 2]; // espectro do sinal filtrado, em dB

// ---------- Etapa 1: gerar o sinal de teste (3 tons somados) ----------
void gerarSinal(float *saida) {
    static float tone1[N];
    static float tone2[N];
    static float tone3[N];

    dsps_tone_gen_f32(tone1, N, 1.0f, toFreqNorm(f1), 0.0f);
    dsps_tone_gen_f32(tone2, N, 1.0f, toFreqNorm(f2), 0.0f);
    dsps_tone_gen_f32(tone3, N, 1.0f, toFreqNorm(f3), 0.0f);

    dsps_add_f32(tone1, tone2, saida, N, 1, 1, 1);
    dsps_add_f32(saida, tone3, saida, N, 1, 1, 1);
}

// ---------- Etapa 2: filtro passa-faixa ----------
static float coeffs1[5] = {
    0.05150721440245413f, 0.0f, -0.05150721440245413f,
    -1.534693565180896f, 0.8969855711950917f
};
static float coeffs2[5] = {
    0.021998737686764813f, 0.0f, -0.021998737686764813f,
    -1.5824392834631165f, 0.9560025246264705f
};

void filtrarPassaFaixa(const float *entrada, float *saida) {
    static float estagio1[N];
    float w1[2] = {0.0f, 0.0f};
    float w2[2] = {0.0f, 0.0f};
    dsps_biquad_f32(entrada, estagio1, N, coeffs1, w1);
    dsps_biquad_f32(estagio1, saida, N, coeffs2, w2);
}

// ---------- Etapa 3: espectro de magnitude em dB ----------
void calcularEspectroDb(const float *sinal, float *magnitudeDb) {
    static float fft_buf[N * 2];
    static float hann[N];

    dsps_wind_hann_f32(hann, N);
    for (int i = 0; i < N; i++) {
        fft_buf[2 * i + 0] = sinal[i] * hann[i];
        fft_buf[2 * i + 1] = 0.0f;
    }

    dsps_fft2r_fc32(fft_buf, N);
    dsps_bit_rev_fc32(fft_buf, N);
    dsps_cplx2reC_fc32(fft_buf, N);

    for (int i = 0; i < N / 2; i++) {
        const float re = fft_buf[2 * i + 0];
        const float im = fft_buf[2 * i + 1];
        magnitudeDb[i] = 10 * log10f((re * re + im * im) / (N * N) + 1e-12f);
    }
}

// ---------- Etapa 4: enviar tudo para visualizacao ----------
void plotarResultados(const float *sinal, const float *filtrado,
                      const float *magnitudeDb, const float *magnitudeFiltradoDb) {
    wserial.println("Espectro (FFT) do sinal 100 Hz + 800 Hz + 2.2 kHz @ fs=8 kHz");

    wserial.plot("signal", 1, sinal, N, "V");
    wserial.plot("filtered", 1, filtrado, N, "V");
    wserial.plot("magnitude", 1, magnitudeDb, N / 2, "dB");
    wserial.plot("magnitude_filtered", 1, magnitudeFiltradoDb, N / 2, "dB");

    dsps_view(magnitudeDb, N / 2, 128, 20, -80, 0, '*');
}

void setup() {
    wserial.begin(115200);
    delay(1000);

    esp_err_t ret = dsps_fft2r_init_fc32(NULL, N);
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

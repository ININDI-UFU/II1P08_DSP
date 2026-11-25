// =============================================================================
//  main0.cpp - Filtro digital passa-faixa + FFT no ESP32 (biblioteca esp-dsp)
// -----------------------------------------------------------------------------
//  Roteiro (executado uma unica vez, no setup):
//    1. Gera um sinal de teste x[n] = soma de 3 senos (100 Hz + 800 Hz + 2200 Hz)
//    2. Filtra x[n] com um passa-faixa IIR (biquads em cascata) -> y[n]
//       O filtro deixa passar 800 Hz e atenua 100 Hz e 2200 Hz.
//    3. Calcula o espectro (FFT, em dB) de x[n] e de y[n]
//    4. Envia os 4 graficos pela serial (LasecPlot / Teleplot)
//  Os coeficientes do filtro sao gerados pelo script Python
//  python/findDigitalFilterIIRBiquadsByTargetFreq.py (veja o README.md).
// =============================================================================
#include <Arduino.h>
#include "services/wserial.h"  // envio de graficos pela serial
#include "dsps_tone_gen.h"     // esp-dsp: gerador de senos (tons)
#include "dsps_add.h"          // esp-dsp: soma de vetores
#include "dsps_biquad.h"       // esp-dsp: filtro IIR de 2a ordem (biquad)
#include "dsps_fft2r.h"        // esp-dsp: FFT radix-2
#include "dsps_wind_hann.h"    // esp-dsp: janela de Hann

// ---------- Parametros do sinal ----------
constexpr int   N  = 1024;     // numero de amostras (a FFT exige potencia de 2)
constexpr float fs = 8000.0f;  // frequencia de amostragem [Hz]
constexpr float f1 = 100.0f;   // tom abaixo da banda -> deve ser atenuado
constexpr float f2 = 800.0f;   // tom dentro da banda -> deve passar
constexpr float f3 = 2200.0f;  // tom acima da banda  -> deve ser atenuado

// ---------- Coeficientes do filtro (copiados da saida do Python) ----------
// Python: fDesejada=[700, 900], ordem=2, fs=8000, filterType="bandpass", desvio=0.05
// Cada linha e um biquad {b0, b1, b2, a1, a2} (a0 = 1), que calcula:
//   y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

float biquads_df2t[NUM_BIQUADS][5] = {
    {1.553585316e-02f, 3.107170631e-02f, 1.553585316e-02f, -1.345224238e+00f, 8.032745419e-01f},
    {1.000000000e+00f, -2.000000000e+00f, 1.000000000e+00f, -1.590652875e+00f, 8.464108338e-01f},
};
constexpr int NUM_BIQUADS = sizeof(biquads_df2t) / sizeof(biquads_df2t[0]);  // = numero de linhas

// ---------- Vetores de amostras ----------
float x[N];         // sinal original
float y[N];         // sinal filtrado
float janela[N];    // janela de Hann (aplicada antes da FFT)
float X_dB[N / 2];  // espectro de x em dB (raias de 0 ate fs/2)
float Y_dB[N / 2];  // espectro de y em dB

// 1. Gera sinal[n] = sen(2*pi*f1*n/fs) + sen(2*pi*f2*n/fs) + sen(2*pi*f3*n/fs)
//    dsps_tone_gen_f32(saida, N, amplitude, f/fs, fase em graus)
//    Atencao: a frequencia e normalizada por fs (f/fs), e NAO por fs/2,
//    apesar de a documentacao da esp-dsp falar em "Nyquist -1..1".
void gerarSinal(float *sinal) {
    static float tom[N];                            // vetor auxiliar: um tom por vez
    dsps_tone_gen_f32(sinal, N, 1.0f, f1 / fs, 0);  // sinal = tom de f1
    dsps_tone_gen_f32(tom, N, 1.0f, f2 / fs, 0);    // tom   = tom de f2
    dsps_add_f32(sinal, tom, sinal, N, 1, 1, 1);    // sinal = sinal + tom
    dsps_tone_gen_f32(tom, N, 1.0f, f3 / fs, 0);    // tom   = tom de f3
    dsps_add_f32(sinal, tom, sinal, N, 1, 1, 1);    // sinal = sinal + tom
}

// 2. Passa o sinal por todos os biquads, um depois do outro (cascata)
void filtrar(const float *entrada, float *saida) {
    memcpy(saida, entrada, N * sizeof(float));  // saida comeca como copia da entrada
    for (int k = 0; k < NUM_BIQUADS; k++) {
        float w[2] = {0, 0};                                // memoria do biquad k (comeca zerada)
        dsps_biquad_f32(saida, saida, N, biquads_df2t[k], w);  // saida = biquad_k(saida)
    }
}

// 3. Espectro em dB: janela -> FFT -> modulo |X[k]| -> dB
void calcularEspectroDb(const float *sinal, float *magDb) {
    static float buf[2 * N];  // complexos intercalados: re0, im0, re1, im1, ...
    for (int n = 0; n < N; n++) {
        buf[2 * n]     = sinal[n] * janela[n];  // parte real = sinal janelado
        buf[2 * n + 1] = 0;                     // parte imaginaria = 0 (sinal real)
    }
    dsps_fft2r_fc32(buf, N);    // FFT (as raias saem embaralhadas, em "bit reverso")
    dsps_bit_rev_fc32(buf, N);  // reordena as raias: 0, 1, 2, ..., N-1
    for (int k = 0; k < N / 2; k++) {  // so ate fs/2: a outra metade e espelhada
        float re = buf[2 * k], im = buf[2 * k + 1];
        magDb[k] = 20 * log10f(sqrtf(re * re + im * im) / N + 1e-6f);  // +1e-6 evita log(0)
    }
}

void setup() {
    wserial.begin(921600);  // mesma velocidade do monitor_speed do platformio.ini
    delay(1000);

    dsps_fft2r_init_fc32(NULL, N);  // prepara a tabela de senos/cossenos da FFT
    dsps_wind_hann_f32(janela, N);  // calcula a janela de Hann (uma unica vez)

    gerarSinal(x);                // 1. x[n]
    filtrar(x, y);                // 2. y[n]
    calcularEspectroDb(x, X_dB);  // 3. |X| em dB
    calcularEspectroDb(y, Y_dB);  //    |Y| em dB

    // 4. Graficos. Eixo horizontal = indice: amostra n (tempo) ou raia k (f = k*fs/N)
    wserial.plot("x", 1, x, N, "V");
    wserial.plot("y", 1, y, N, "V");
    wserial.plot("X_dB", 1, X_dB, N / 2, "dB");
    wserial.plot("Y_dB", 1, Y_dB, N / 2, "dB");
}

void loop() {
    wserial.update();  // nada a processar: so mantem a serial atendida
}

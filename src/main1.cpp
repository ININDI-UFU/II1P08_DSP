#include <Arduino.h>
#include <math.h>

#include "services/wserial.h"
#include "dsps_biquad.h"

// Versao simulada (sem hardware) do main2.cpp:
// - gera um sinal periodico sintetico no lugar da leitura real do ADC
// - aplica os mesmos dois biquads passa-baixa em cascata
// - envia sinal original e sinal filtrado pela serial (wserial/Teleplot)

constexpr uint32_t SERIAL_BAUD = 921600;
constexpr float SAMPLE_RATE_HZ = 100.0f;
constexpr uint32_t SAMPLE_PERIOD_US = static_cast<uint32_t>(1000000.0f / SAMPLE_RATE_HZ);

// Sinal simulado: 10 Hz (deve passar) + 40 Hz (deve ser atenuado pelo filtro abaixo).
// As duas fecham exatamente em 10 amostras @ fs=100 Hz:
// 10 Hz -> 1 ciclo em 10 amostras | 40 Hz -> 4 ciclos em 10 amostras.
constexpr float SIGNAL_FREQ_PASS = 10.0f;
constexpr float SIGNAL_FREQ_STOP = 40.0f;
constexpr int   SIGNAL_PERIOD    = 10;

static float signal_lut[SIGNAL_PERIOD];
static int lut_index = 0;

// Butterworth passa-baixa de 4a ordem, fs=100 Hz, fc=20 Hz.
// Esta frequencia de corte deixa 10 Hz passar bem e atenua fortemente 40 Hz.
// Coeficientes no formato exigido por dsps_biquad_f32:
// [b0, b1, b2, a1, a2], com a0 = 1.
static float lowpass_stage1_coeffs[5] = {
    0.046582906636443676f,
    0.09316581327288735f,
    0.046582906636443676f,
    -0.3289756779965208f,
    0.06458765330651409f,
};

static float lowpass_stage2_coeffs[5] = {
    1.0f,
    2.0f,
    1.0f,
    -0.4531195201603722f,
    0.4663255708081653f,
};

static float lowpass_stage1_state[2] = {0.0f, 0.0f};
static float lowpass_stage2_state[2] = {0.0f, 0.0f};

void prepararSinalSimulado() {
    for (int i = 0; i < SIGNAL_PERIOD; i++) {
        const float n = static_cast<float>(i);
        signal_lut[i] =
            sinf(2.0f * PI * SIGNAL_FREQ_PASS * n / SAMPLE_RATE_HZ) +
            sinf(2.0f * PI * SIGNAL_FREQ_STOP * n / SAMPLE_RATE_HZ);
    }
}

float readSimulatedVolts() {
    const float sample = signal_lut[lut_index];
    lut_index++;
    if (lut_index >= SIGNAL_PERIOD) {
        lut_index = 0;
    }
    return sample;
}

float filterLowpassCascade(float sample) {
    float stage1 = 0.0f;
    float filtered = 0.0f;

    dsps_biquad_f32(&sample, &stage1, 1, lowpass_stage1_coeffs, lowpass_stage1_state);
    dsps_biquad_f32(&stage1, &filtered, 1, lowpass_stage2_coeffs, lowpass_stage2_state);

    return filtered;
}

void publishSample(float raw, float filtered) {
    wserial.plot("raw", raw, "V");
    wserial.plot("filtered", filtered, "V");
}

void updateAll() {
    const float raw = readSimulatedVolts();
    const float filtered = filterLowpassCascade(raw);
    publishSample(raw, filtered);
}

void setup() {
    wserial.begin(SERIAL_BAUD);
    delay(1000);
    prepararSinalSimulado();
    wserial.println("Simulacao @ 100 Hz: sinal sintetico e passa-baixa biquad em cascata");
}

void loop() {
    static uint32_t lastMs = 0;
    const uint32_t now = micros();

    if ((now - lastMs) >= SAMPLE_PERIOD_US) {
        lastMs = now;
        updateAll(); // ponto de entrada principal para breakpoints
    }

    wserial.update();
}

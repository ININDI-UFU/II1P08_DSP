#include <Arduino.h>

#include "services/wserial.h"
#include "dsps_biquad.h"

// Aquisicao real no GPIO34 @ 100 Hz:
// - le o ADC uma vez a cada 10 ms
// - aplica dois biquads passa-baixa em cascata
// - envia sinal original e sinal filtrado pela serial (wserial/Teleplot)

constexpr uint8_t ADC_PIN = 34;
constexpr uint32_t SERIAL_BAUD = 921600;

constexpr float SAMPLE_RATE_HZ = 100.0f;
constexpr uint32_t SAMPLE_PERIOD_US = static_cast<uint32_t>(1000000.0f / SAMPLE_RATE_HZ);
constexpr float ADC_MAX_COUNTS = 4095.0f;
constexpr float ADC_REF_VOLTS = 3.3f;

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

float readAdcVolts() {
    const uint16_t counts = analogRead(ADC_PIN);
    return (static_cast<float>(counts) * ADC_REF_VOLTS) / ADC_MAX_COUNTS;
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
    const float raw = readAdcVolts();
    const float filtered = filterLowpassCascade(raw);
    publishSample(raw, filtered);
}

void setup() {
    wserial.begin(SERIAL_BAUD);
    delay(1000);
    analogReadResolution(12);
    analogSetPinAttenuation(ADC_PIN, ADC_11db);
    pinMode(ADC_PIN, INPUT);
    wserial.println("GPIO34 @ 100 Hz: sinal original e passa-baixa biquad em cascata");
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

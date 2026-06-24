// Leitura de um sinal analógico via ADC contínuo (DMA) e cálculo do TRUE RMS
// por janela, em tempo real, usando esp-dsp.
//
// A biblioteca esp-dsp não possui uma função dsps_rms (nem a versão empacotada
// com este framework, nem a atual em github.com/espressif/esp-dsp). O true RMS
// é obtido pela própria definição matemática: RMS = sqrt(mean(x[i]^2)), onde a
// soma dos quadrados é acelerada em hardware por dsps_dotprod_f32
// (dest = soma(src1[i]*src2[i])), passando o sinal como os dois operandos.
#include <Arduino.h>
#include <math.h>

#include "services/wserial.h"
#include "services/AdcDmaEsp.h"
#include "dsps_dotprod.h"

// ---------- Parâmetros ----------
constexpr adc_channel_t ADC_CH      = ADC_CHANNEL_3;   // GPIO4 no ESP32-S3
constexpr uint32_t      SAMPLE_HZ   = 20000;           // taxa de amostragem do ADC [Hz]
constexpr int           WINDOW_LEN  = 256;             // amostras por janela de RMS (~12.8 ms @ 20 kHz)
constexpr float         ADC_TO_VOLT = 3.3f / 4095.0f;  // 12 bits, ATTEN_DB_12 -> faixa 0-3.3V

static float window[WINDOW_LEN];
static int   windowIdx = 0;

void setup() {
    wserial.begin(115200);
    delay(1000);

    if (!adcDma.begin(ADC_CH, SAMPLE_HZ)) {
        wserial.println("Falha ao inicializar o ADC contínuo (DMA)");
        return;
    }
    adcDma.start();
    wserial.println("RMS em tempo real iniciado");
}

void loop() {
    wserial.update();

    if (!adcDma.available()) return;

    AdcDmaSample buf[ADC_DMA_SAMPLES_PER_FRAME];
    const size_t n = adcDma.read(buf, ADC_DMA_SAMPLES_PER_FRAME);

    for (size_t i = 0; i < n; i++) {
        window[windowIdx++] = buf[i].value * ADC_TO_VOLT;

        if (windowIdx == WINDOW_LEN) {
            float sumSq = 0.0f;
            dsps_dotprod_f32(window, window, &sumSq, WINDOW_LEN); // sumSq = soma(window[i]^2)
            const float rms = sqrtf(sumSq / WINDOW_LEN);

            wserial.plot("rms", rms, "V");

            windowIdx = 0; // próxima janela (não sobreposta)
        }
    }
}

// Mesmo sinal do main0.cpp (soma de 3 senoides: 100 Hz + 800 Hz + 2.2 kHz,
// fs=8 kHz, filtro passa-faixa isolando o tom de 800 Hz) só que gerado de
// forma CONTÍNUA: aqui não existe buffer de N amostras computado uma vez em
// setup() e plotado de uma vez só. Em vez disso, cada amostra é sintetizada
// e filtrada individualmente dentro de loop(), no instante certo (acumulador
// de fase + filtro IIR com estado persistente entre chamadas), exatamente
// como um gerador de função real geraria a forma de onda continuamente --
// só que, em vez de sair num pino DAC, ela é escrita amostra a amostra na
// serial, para sempre.
#include <Arduino.h>
#include <math.h>

#include "services/wserial.h"
#include "dsps_biquad.h"
#include "dsps_biquad_gen.h"

// ---------- Parâmetros do sinal (iguais ao main0.cpp) ----------
constexpr float fs = 8000.0f; // frequência de amostragem [Hz]
constexpr float f1 = 100.0f;  // tom 1 [Hz] (mais distante da banda de 800 Hz)
constexpr float f2 = 800.0f;  // tom 2 [Hz] -> é o tom que o filtro passa-faixa isola
constexpr float f3 = 2200.0f; // tom 3 [Hz] (mais distante da banda de 800 Hz)

constexpr uint32_t SAMPLE_PERIOD_US = (uint32_t)(1000000.0f / fs); // 125 us @ 8 kHz

// A serial a 115200 bps aguenta uns poucos milhares de linhas/s, não as 8000
// amostras/s geradas internamente -- por isso só 1 a cada DECIMATION amostras
// é enviada para a serial (a geração e o filtro continuam correndo a fs=8 kHz
// "por debaixo", só a publicação é que é mais lenta). Resultado: ~100 Hz de
// taxa de atualização no plot, suficiente para acompanhar visualmente os 3
// tons (o maior, 1.5 kHz, ainda fica bem caracterizado por amostragem).
constexpr uint32_t DECIMATION = fs / 100; // publica a ~100 amostras/s

// Incremento de fase por amostra (rad/amostra) de cada tom -- equivalente ao
// "freq_norm" do main0.cpp, só que aqui aplicado amostra a amostra.
constexpr float dphi1 = 2.0f * PI * f1 / fs;
constexpr float dphi2 = 2.0f * PI * f2 / fs;
constexpr float dphi3 = 2.0f * PI * f3 / fs;

// ---------- Estado persistente entre chamadas de loop() ----------
static float phase1 = 0.0f, phase2 = 0.0f, phase3 = 0.0f; // fase corrente de cada tom
static float bpfCoeffs[5];                                // coeficientes do biquad (calculados 1x em setup)
static float w1[2] = {0.0f, 0.0f};                         // estado do 1º estágio do filtro
static float w2[2] = {0.0f, 0.0f};                         // estado do 2º estágio (em cascata)
static uint32_t nextSampleUs = 0;                          // instante (us) da próxima amostra
static uint32_t sampleCount  = 0;                          // conta amostras geradas, p/ decimação

// ---------- Gera 1 amostra do sinal somado, avançando as 3 fases ----------
float gerarAmostra() {
    const float amostra = sinf(phase1) + sinf(phase2) + sinf(phase3);

    phase1 += dphi1; if (phase1 >= 2.0f * PI) phase1 -= 2.0f * PI;
    phase2 += dphi2; if (phase2 >= 2.0f * PI) phase2 -= 2.0f * PI;
    phase3 += dphi3; if (phase3 >= 2.0f * PI) phase3 -= 2.0f * PI;

    return amostra;
}

// ---------- Filtra 1 amostra (cascata de 2 biquads, mesmo desenho do main0.cpp) ----------
// dsps_biquad_f32 com len=1 processa uma única amostra por chamada; como w1/w2
// são estáticos (não são zerados a cada chamada), o estado do filtro IIR é
// preservado de uma amostra para a próxima -- é isso que torna o filtro
// "contínuo" em vez de recalcular tudo do zero a cada N amostras.
float filtrarAmostra(float entrada) {
    float estagio1, saida;
    dsps_biquad_f32(&entrada, &estagio1, 1, bpfCoeffs, w1); // 1º estágio
    dsps_biquad_f32(&estagio1, &saida, 1, bpfCoeffs, w2);   // 2º estágio (em série)
    return saida;
}

void setup() {
    wserial.begin(115200);
    delay(1000);

    constexpr float bpfQ = 5.0f;
    dsps_biquad_gen_bpf0db_f32(bpfCoeffs, f2 / fs, bpfQ);

    wserial.println("Gerador contínuo: 100 Hz + 800 Hz + 2.2 kHz @ fs=8 kHz (saída ao vivo na serial)");
    nextSampleUs = micros();
}

void loop() {
    wserial.update();

    // Pacing não bloqueante: só gera a próxima amostra quando o tempo dela
    // chegar (comparação com cast para int32_t é segura mesmo quando micros()
    // dá overflow, pois o resultado do subtração continua correto em wraparound).
    const uint32_t now = micros();
    if ((int32_t)(now - nextSampleUs) < 0) return;
    nextSampleUs += SAMPLE_PERIOD_US;

    const float sinal    = gerarAmostra();
    const float filtrado = filtrarAmostra(sinal);

    if (++sampleCount % DECIMATION == 0) {
        wserial.plot("signal", sinal, "V");
        wserial.plot("filtered", filtrado, "V");
    }
}

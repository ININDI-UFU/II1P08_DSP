#include <Arduino.h>
#include <math.h>

#include "services/wserial.h"
#include "dsps_fft2r.h"
#include "dsps_view.h"
#include "dsps_wind_hann.h"

// Simulacao de aquisicao continua tipo ADC+DMA:
// - um gerador periodico alimenta frames pequenos ("DMA")
// - os frames entram em buffers circulares
// - o filtro roda continuamente, preservando estado entre amostras
// - a FFT usa snapshots da janela mais recente
// - a serial publica em lotes e em taxa menor para o plot acompanhar

// ---------- Parametros do sinal ----------
constexpr int   N  = 1024;    // janela da FFT
constexpr float fs = 8000.0f; // taxa de amostragem
constexpr float f1 = 100.0f;
constexpr float f2 = 800.0f;
constexpr float f3 = 2200.0f;

// Todas as tres senoides fecham exatamente em 80 amostras:
// 100 Hz  -> 1 ciclo em 80 amostras
// 800 Hz  -> 8 ciclos em 80 amostras
// 2200 Hz -> 22 ciclos em 80 amostras
constexpr int SIGNAL_PERIOD = 80;

// Simula o frame entregue pelo DMA a cada interrupcao.
constexpr int DMA_FRAME_SAMPLES = 64;

// Processa DSP em janelas deslizantes.
constexpr int PROCESS_HOP = 256;

// Publica mais devagar do que processa para nao saturar a serial.
constexpr uint32_t PUBLISH_INTERVAL_MS = 1000;

// Mantem poucos pontos no dominio do tempo para o plot suportar bem.
constexpr int TIME_PLOT_POINTS = 128;
constexpr int FFT_PLOT_POINTS  = 128;
constexpr int FFT_DECIMATION   = (N / 2) / FFT_PLOT_POINTS;
constexpr uint32_t SERIAL_BAUD = 921600;

static_assert((N & (N - 1)) == 0, "N deve ser potencia de 2.");
static_assert((N / 2) % FFT_PLOT_POINTS == 0, "A decimacao da FFT deve ser inteira.");

// ---------- Buffers ----------
static float signal_lut[SIGNAL_PERIOD];
static float raw_ring[N];
static float filtered_ring[N];

static float raw_frame[N];
static float filtered_frame[N];
static float magnitude[N / 2];
static float magnitude_filtered[N / 2];

static float plot_raw[TIME_PLOT_POINTS];
static float plot_filtered[TIME_PLOT_POINTS];
static float plot_magnitude[FFT_PLOT_POINTS];
static float plot_magnitude_filtered[FFT_PLOT_POINTS];

static float hann[N];

// ---------- Estado da simulacao ----------
static int lut_index = 0;
static int ring_write_index = 0;
static uint32_t samples_since_process = 0;
static uint32_t processed_windows = 0;
static bool dsp_snapshot_ready = false;

// Estado do biquad em forma direta II transposta, preservado continuamente.
static float stage1_state[2] = {0.0f, 0.0f};
static float stage2_state[2] = {0.0f, 0.0f};

// Coeficientes [b0, b1, b2, a1, a2], a0 = 1.
static const float coeffs1[5] = {
    0.05150721440245413f, 0.0f, -0.05150721440245413f,
    -1.534693565180896f, 0.8969855711950917f
};
static const float coeffs2[5] = {
    0.021998737686764813f, 0.0f, -0.021998737686764813f,
    -1.5824392834631165f, 0.9560025246264705f
};

void prepararLutPeriodica() {
    for (int i = 0; i < SIGNAL_PERIOD; i++) {
        const float n = static_cast<float>(i);
        signal_lut[i] =
            sinf(2.0f * PI * f1 * n / fs) +
            sinf(2.0f * PI * f2 * n / fs) +
            sinf(2.0f * PI * f3 * n / fs);
    }
}

inline float lerLutPeriodica() {
    const float sample = signal_lut[lut_index];
    lut_index++;
    if (lut_index >= SIGNAL_PERIOD) {
        lut_index = 0;
    }
    return sample;
}

float aplicarBiquadAmostra(float x, const float *c, float state[2]) {
    const float y  = c[0] * x + state[0];
    const float s0 = c[1] * x - c[3] * y + state[1];
    const float s1 = c[2] * x - c[4] * y;
    state[0] = s0;
    state[1] = s1;
    return y;
}

float filtrarAmostraContinua(float sample) {
    const float y1 = aplicarBiquadAmostra(sample, coeffs1, stage1_state);
    return aplicarBiquadAmostra(y1, coeffs2, stage2_state);
}

void armazenarAmostra(float raw_sample) {
    const float filtered_sample = filtrarAmostraContinua(raw_sample);
    raw_ring[ring_write_index] = raw_sample;
    filtered_ring[ring_write_index] = filtered_sample;

    ring_write_index++;
    if (ring_write_index >= N) {
        ring_write_index = 0;
    }
}

void copiarJanelaDoRing(const float *ring, float *frame) {
    for (int i = 0; i < N; i++) {
        const int idx = ring_write_index + i;
        frame[i] = ring[idx < N ? idx : idx - N];
    }
}

void calcularEspectroDb(const float *sinal, float *magnitudeDb) {
    static float fft_buf[N * 2];

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
        magnitudeDb[i] = 10.0f * log10f((re * re + im * im) / (N * N) + 1e-12f);
    }
}

void atualizarSnapshotDsp() {
    copiarJanelaDoRing(raw_ring, raw_frame);
    copiarJanelaDoRing(filtered_ring, filtered_frame);
    calcularEspectroDb(raw_frame, magnitude);
    calcularEspectroDb(filtered_frame, magnitude_filtered);
    dsp_snapshot_ready = true;
    processed_windows++;
}

void processarSeJanelaPronta() {
    while (samples_since_process >= PROCESS_HOP) {
        samples_since_process -= PROCESS_HOP;
        atualizarSnapshotDsp();
    }
}

void simularFrameDma() {
    for (int i = 0; i < DMA_FRAME_SAMPLES; i++) {
        armazenarAmostra(lerLutPeriodica());
    }
    samples_since_process += DMA_FRAME_SAMPLES;
    processarSeJanelaPronta();
}

void atualizarAquisicaoContinua() {
    static uint32_t last_us = micros();
    static uint32_t sample_accumulator_us = 0;
    static uint32_t pending_samples = 0;

    const uint32_t now_us = micros();
    const uint32_t elapsed_us = now_us - last_us;
    last_us = now_us;
    sample_accumulator_us += elapsed_us;

    constexpr uint32_t SAMPLE_PERIOD_US = static_cast<uint32_t>(1000000.0f / fs);
    pending_samples += sample_accumulator_us / SAMPLE_PERIOD_US;
    sample_accumulator_us %= SAMPLE_PERIOD_US;

    while (pending_samples >= DMA_FRAME_SAMPLES) {
        simularFrameDma();
        pending_samples -= DMA_FRAME_SAMPLES;
    }
}

void prepararCurvasParaPlot() {
    constexpr int TIME_PLOT_START = N - TIME_PLOT_POINTS;

    for (int i = 0; i < TIME_PLOT_POINTS; i++) {
        const int src = TIME_PLOT_START + i;
        plot_raw[i] = raw_frame[src];
        plot_filtered[i] = filtered_frame[src];
    }

    for (int i = 0; i < FFT_PLOT_POINTS; i++) {
        const int src = i * FFT_DECIMATION;
        plot_magnitude[i] = magnitude[src];
        plot_magnitude_filtered[i] = magnitude_filtered[src];
    }
}

void publicarSnapshot() {
    if (!dsp_snapshot_ready) {
        return;
    }

    prepararCurvasParaPlot();

    wserial.print("#batch:");
    wserial.println(processed_windows);
    wserial.plot("signal", 1, plot_raw, TIME_PLOT_POINTS, "V");
    wserial.plot("filtered", 1, plot_filtered, TIME_PLOT_POINTS, "V");
    wserial.plot("magnitude", 1, plot_magnitude, FFT_PLOT_POINTS, "dB");
    wserial.plot("magnitude_filtered", 1, plot_magnitude_filtered, FFT_PLOT_POINTS, "dB");
}

void iniciarPreenchimentoDoRing() {
    for (int i = 0; i < N; i++) {
        armazenarAmostra(lerLutPeriodica());
    }
    atualizarSnapshotDsp();
}

void setup() {
    wserial.begin(SERIAL_BAUD);
    delay(1000);

    const esp_err_t ret = dsps_fft2r_init_fc32(NULL, N);
    if (ret != ESP_OK) {
        wserial.println(String("Falha ao inicializar FFT: ") + ret);
        return;
    }

    dsps_wind_hann_f32(hann, N);
    prepararLutPeriodica();
    iniciarPreenchimentoDoRing();

    wserial.println("Simulacao continua ADC+DMA: filtro e FFT por janela deslizante");
}

void loop() {
    atualizarAquisicaoContinua();

    static uint32_t last_publish_ms = 0;
    if (millis() - last_publish_ms >= PUBLISH_INTERVAL_MS) {
        last_publish_ms = millis();
        publicarSnapshot();
    }

    wserial.update();
}

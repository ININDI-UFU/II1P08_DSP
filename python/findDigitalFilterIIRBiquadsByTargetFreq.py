import argparse
from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, sosfreqz


@dataclass(frozen=True)
class IIRBiquadResult:
    fDesejada: float
    ordem: int
    fs: float
    filterType: str
    desvio: float
    isBP: bool
    fc: float
    f_complementar_hz: float
    magnitude_desejada: float
    magnitude_obtida: float
    erro_magnitude: float
    sos: np.ndarray

    @property
    def quantidade_biquads(self):
        return self.ordem // 2

    @property
    def biquads_df2t(self):
        """Coeficientes [b0, b1, b2, a1, a2] para Direct Form II Transposto."""
        return self.sos[:, [0, 1, 2, 4, 5]]


def findDigitalFilterIIRBiquadsByTargetFreq(
    fDesejada,
    ordem,
    fs,
    filterType="lowpass",
    desvio=0.05,
    isBP=True,
    plot=True,
    verbose=True,
    worN=4096,
):
    """
    Calcula um filtro IIR Butterworth e retorna suas secoes biquad.

    A ordem deve ser par e maior ou igual a 2. O numero de biquads retornado
    sera ordem / 2.

    Retorno principal:
      resultado.sos          -> [b0, b1, b2, a0, a1, a2]
      resultado.biquads_df2t -> [b0, b1, b2, a1, a2]

    isBP=True significa que fDesejada esta na banda de passagem.
    isBP=False significa que fDesejada esta na banda de rejeicao.
    """
    fDesejada, ordem, fs, filterType, desvio = _validar_parametros(
        fDesejada=fDesejada,
        ordem=ordem,
        fs=fs,
        filterType=filterType,
        desvio=desvio,
    )

    magnitude_desejada = 1 - desvio if isBP else desvio
    fc = _calcular_fc_digital(
        fDesejada=fDesejada,
        ordem=ordem,
        fs=fs,
        filterType=filterType,
        magnitude=magnitude_desejada,
    )

    magnitude_complementar = desvio if isBP else 1 - desvio
    f_complementar_hz = _calcular_frequencia_por_magnitude(
        fc=fc,
        ordem=ordem,
        fs=fs,
        filterType=filterType,
        magnitude=magnitude_complementar,
    )

    sos = butter(ordem, fc, btype=filterType, fs=fs, output="sos")
    magnitude_obtida = _magnitude_exata(sos, fDesejada, fs)

    resultado = IIRBiquadResult(
        fDesejada=fDesejada,
        ordem=ordem,
        fs=fs,
        filterType=filterType,
        desvio=desvio,
        isBP=bool(isBP),
        fc=fc,
        f_complementar_hz=f_complementar_hz,
        magnitude_desejada=magnitude_desejada,
        magnitude_obtida=magnitude_obtida,
        erro_magnitude=abs(magnitude_obtida - magnitude_desejada),
        sos=sos,
    )

    if verbose:
        _imprimir_resultado(resultado)

    if plot:
        _plotar_resultado(resultado, worN)

    return resultado


def formatBiquadsForC(resultado, nome_array="biquads_df2t", incluir_a0=False, casas=9):
    """
    Formata os biquads em um array C.

    Por padrao usa [b0, b1, b2, a1, a2], ideal para Direct Form II Transposto.
    Use incluir_a0=True para gerar [b0, b1, b2, a0, a1, a2].
    """
    matriz = resultado.sos if incluir_a0 else resultado.biquads_df2t
    linhas = [f"const float {nome_array}[{matriz.shape[0]}][{matriz.shape[1]}] = {{"]

    for secao in matriz:
        valores = ", ".join(f"{valor:.{casas}e}f" for valor in secao)
        linhas.append(f"    {{{valores}}},")

    linhas.append("};")
    return "\n".join(linhas)


def _validar_parametros(fDesejada, ordem, fs, filterType, desvio):
    fDesejada = float(fDesejada)
    ordem = int(ordem)
    fs = float(fs)
    filterType = str(filterType).lower()
    desvio = float(desvio)

    if filterType not in ("lowpass", "highpass"):
        raise ValueError("filterType deve ser 'lowpass' ou 'highpass'.")

    if fs <= 0:
        raise ValueError("fs deve ser maior que zero.")

    if not (0 < fDesejada < fs / 2):
        raise ValueError(f"fDesejada={fDesejada} Hz deve estar entre 0 e fs/2 ({fs / 2} Hz).")

    if ordem < 2 or ordem % 2 != 0:
        raise ValueError("ordem deve ser par e maior ou igual a 2.")

    if not (0 < desvio < 1):
        raise ValueError("desvio deve estar entre 0 e 1.")

    return fDesejada, ordem, fs, filterType, desvio


def _calcular_fator_butterworth(magnitude, ordem):
    if not (0 < magnitude < 1):
        raise ValueError("A magnitude deve estar entre 0 e 1.")

    return ((1 / (magnitude**2)) - 1) ** (1 / (2 * ordem))


def _calcular_fc_digital(fDesejada, ordem, fs, filterType, magnitude):
    fator = _calcular_fator_butterworth(magnitude, ordem)
    wd = 2 * np.pi * fDesejada / fs
    tan_wd = np.tan(wd / 2)

    if filterType == "lowpass":
        tan_wc = tan_wd / fator
    else:
        tan_wc = tan_wd * fator

    wc = 2 * np.arctan(tan_wc)
    fc = wc * fs / (2 * np.pi)
    nyq = fs / 2

    if not (0 < fc < nyq):
        raise ValueError(f"fc={fc:.6g} Hz ficou fora da faixa valida (0, {nyq:.6g} Hz).")

    return float(fc)


def _calcular_frequencia_por_magnitude(fc, ordem, fs, filterType, magnitude):
    fator = _calcular_fator_butterworth(magnitude, ordem)
    wc = 2 * np.pi * fc / fs
    tan_wc = np.tan(wc / 2)

    if filterType == "lowpass":
        tan_w = tan_wc * fator
    else:
        tan_w = tan_wc / fator

    w = 2 * np.arctan(tan_w)
    return float(w * fs / (2 * np.pi))


def _magnitude_exata(sos, frequencia_hz, fs):
    w = [2 * np.pi * frequencia_hz / fs]
    _, h = sosfreqz(sos, worN=w)
    return float(np.abs(h[0]))


def _imprimir_resultado(resultado):
    print("Resultado do filtro IIR em biquads:")
    print(f"  Tipo: {resultado.filterType}")
    print(f"  Ordem: {resultado.ordem}")
    print(f"  Quantidade de biquads: {resultado.quantidade_biquads}")
    print(f"  Frequencia desejada: {resultado.fDesejada:.6f} Hz")
    print(f"  Frequencia de corte: {resultado.fc:.6f} Hz")
    print(f"  Frequencia complementar: {resultado.f_complementar_hz:.6f} Hz")
    print(f"  Magnitude alvo/obtida: {resultado.magnitude_desejada:.9f} / {resultado.magnitude_obtida:.9f}")


def _imprimir_sos(resultado):
    print("  SOS [b0, b1, b2, a0, a1, a2]:")
    for i, secao in enumerate(resultado.sos, start=1):
        valores = ", ".join(f"{coef:.9e}" for coef in secao)
        print(f"    biquad_{i}: {valores}")


def _plotar_resultado(resultado, worN):
    from matplotlib.pyplot import show, subplots

    w, h = sosfreqz(resultado.sos, worN=worN)
    freqs_hz = w * resultado.fs / (2 * np.pi)
    h_mag = np.abs(h)

    fig, ax = subplots(figsize=(14, 5))
    ax.plot(freqs_hz, h_mag, label="Resposta do filtro IIR em biquads", color="blue")
    ax.axvline(resultado.fDesejada, color="yellow", linestyle="--", label="Frequencia desejada")
    ax.axvline(resultado.fc, color="green", linestyle="--", label="Frequencia de corte")
    ax.axvline(resultado.f_complementar_hz, color="purple", linestyle="--", label="Frequencia complementar")
    ax.scatter([resultado.fDesejada], [resultado.magnitude_desejada], color="red", label="Magnitude alvo")

    f_min_plot = max(0, min(resultado.fDesejada, resultado.fc, resultado.f_complementar_hz) * 0.8)
    f_max_plot = min(
        resultado.fs / 2,
        max(resultado.fDesejada, resultado.fc, resultado.f_complementar_hz) * 1.5,
    )
    ax.set_xlim(f_min_plot, f_max_plot)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("Frequencia (Hz)")
    ax.set_ylabel("Magnitude")
    ax.set_title(f"Filtro IIR Butterworth - {resultado.quantidade_biquads} biquad(s)")
    ax.grid(True)
    ax.legend()
    show()


def _criar_parser():
    parser = argparse.ArgumentParser(
        description="Calcula filtro IIR Butterworth em biquads a partir de uma frequencia alvo.",
    )
    parser.add_argument("--fDesejada", type=float, default=100.0, help="Frequencia desejada em Hz.")
    parser.add_argument("--ordem", type=int, default=10, help="Ordem do filtro. Deve ser par.")
    parser.add_argument("--fs", type=float, default=1600.0, help="Frequencia de amostragem em Hz.")
    parser.add_argument(
        "--filterType",
        choices=("lowpass", "highpass"),
        default="lowpass",
        help="Tipo do filtro.",
    )
    parser.add_argument("--desvio", type=float, default=0.05, help="Desvio permitido em magnitude.")
    parser.add_argument(
        "--isBP",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="True se fDesejada esta na banda de passagem; False se esta na banda de rejeicao.",
    )
    parser.add_argument(
        "--plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exibe o grafico da resposta em frequencia.",
    )
    parser.add_argument(
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Exibe o resumo do resultado no console.",
    )
    parser.add_argument(
        "--show-sos",
        action="store_true",
        help="Exibe todas as secoes SOS no console.",
    )
    parser.add_argument(
        "--show-c",
        action="store_true",
        help="Exibe os coeficientes formatados para C no console.",
    )
    parser.add_argument(
        "--c-only",
        action="store_true",
        help="Exibe somente os coeficientes em formato C, sem resumo adicional.",
    )
    parser.add_argument(
        "--incluir-a0",
        action="store_true",
        help="Ao usar --show-c, inclui a0 no array gerado.",
    )
    parser.add_argument("--casas", type=int, default=9, help="Casas decimais em notacao cientifica para --show-c.")
    return parser


if __name__ == "__main__":
    args = _criar_parser().parse_args()
    c_only = args.c_only or (not args.verbose and not args.show_sos and not args.show_c)
    verbose = args.verbose and not c_only

    resultado = findDigitalFilterIIRBiquadsByTargetFreq(
        fDesejada=args.fDesejada,
        ordem=args.ordem,
        fs=args.fs,
        filterType=args.filterType,
        desvio=args.desvio,
        isBP=args.isBP,
        plot=args.plot,
        verbose=verbose,
    )

    if args.show_sos and not c_only:
        _imprimir_sos(resultado)

    if args.show_c or c_only:
        if not c_only:
            print("\nCoeficientes para Direct Form II Transposto:")
        print(
            formatBiquadsForC(
                resultado,
                incluir_a0=args.incluir_a0,
                casas=args.casas,
            )
        )

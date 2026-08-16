# -*- coding: utf-8 -*-
"""
Incerteza numérica e portões de validade.

Um número sem barra de erro não é resultado. E um número dentro da barra de erro
mas FORA do regime de validade da teoria com que você compara é pior ainda:
parece certo e está errado.

Este módulo faz as duas coisas:

  1. INCERTEZA NUMÉRICA — extrapolação de Richardson sobre o passo da grade,
     mais o teste de truncamento (o resultado depende de onde eu cortei?).

  2. PORTÕES DE VALIDADE — antes de comparar com qualquer teoria universal,
     conferir se as hipóteses dela valem. |a| >> r0? k*R << 1? r0 << l_ho?

O portão nasceu de um erro real: comparamos o solver confinado com a fórmula de
Busch (alcance zero) em pontos onde |a|/r0 valia 0,6. A fórmula não vale ali, e
a discordância que parecia física era só hipótese violada.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import math


# =============================================================================
#  Medida = valor + incerteza + de onde ela veio + o que pode invalidar
# =============================================================================
@dataclass
class Medida:
    valor: float
    incerteza: float = 0.0
    unidade: str = ""
    orcamento: dict = field(default_factory=dict)   # de onde vem cada pedaço
    alertas: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.alertas

    @property
    def digitos_significativos(self) -> int:
        """Quantos dígitos você tem direito de escrever."""
        if self.incerteza <= 0 or self.valor == 0:
            return 15
        return max(1, int(math.floor(math.log10(abs(self.valor / self.incerteza)))) + 1)

    def __str__(self) -> str:
        if self.incerteza > 0:
            n = max(0, -int(math.floor(math.log10(self.incerteza))) + 1)
            s = f"{self.valor:.{n}f} ± {self.incerteza:.{n}f}"
        else:
            s = f"{self.valor:.10g}"
        if self.unidade:
            s += f" {self.unidade}"
        if self.alertas:
            s += "   [!] " + "; ".join(self.alertas)
        return s

    def compativel_com(self, ref: float, tol_sigma: float = 3.0) -> bool:
        """O valor de referência cai dentro de N sigmas?"""
        if self.incerteza <= 0:
            return math.isclose(self.valor, ref, rel_tol=1e-9)
        return abs(self.valor - ref) <= tol_sigma * self.incerteza


# =============================================================================
#  1. Incerteza de discretização — extrapolação de Richardson
# =============================================================================
def richardson(func, passos, ordem=2, unidade=""):
    """Extrapola func(dr) para dr -> 0 e devolve o resíduo como incerteza.

    func    : callable, recebe o passo e devolve um float
    passos  : lista de passos, do MAIOR para o MENOR (ex.: [4e-3, 2e-3, 1e-3])
    ordem   : ordem esperada do método (2 = diferença central, 4 = Numerov liso)

    A incerteza é |extrapolado - mais_fino|: honesta e conservadora.
    Se a ordem medida discordar da esperada, entra um alerta — é o sinal de
    que a hipótese do método quebrou (ex.: Numerov em potencial descontínuo).
    """
    vals = [float(func(h)) for h in passos]
    alertas = []

    escala = max(abs(v) for v in vals) or 1.0
    if len(vals) >= 3:
        razao = passos[0] / passos[1]
        d1, d2 = vals[0] - vals[1], vals[1] - vals[2]
        if abs(d2) < 1e-13 * escala:
            # ja convergiu ate a precisao de maquina: a "ordem" e so ruido
            p_medida = float("nan")
            d1 = d2 = 0.0
        elif d2 != 0 and d1 / d2 > 0:
            p_medida = math.log(abs(d1 / d2)) / math.log(razao)
            if abs(p_medida - ordem) > 0.5:
                alertas.append(f"ordem medida {p_medida:.1f} != esperada {ordem} "
                               f"(hipótese do método pode ter quebrado)")
        else:
            p_medida = float("nan")
            alertas.append("convergência não monotônica")
    else:
        p_medida = float("nan")

    # Richardson com a ordem esperada
    r = passos[-2] / passos[-1]
    ext = vals[-1] + (vals[-1] - vals[-2]) / (r**ordem - 1)
    inc = abs(ext - vals[-1])

    # se a incerteza ja e desprezivel, a "ordem medida" e ruido: nao alarme
    if inc < 1e-8 * escala:
        alertas = [a for a in alertas if "ordem medida" not in a
                   and "monot" not in a]

    return Medida(valor=ext, incerteza=inc, unidade=unidade,
                  orcamento={"discretizacao": inc, "ordem_medida": p_medida},
                  alertas=alertas)


# =============================================================================
#  2. Incerteza de truncamento — o resultado depende de onde eu cortei?
# =============================================================================
def truncamento(func, cortes, unidade=""):
    """func(corte) para vários r_max / R. A física NÃO pode depender do corte.

    Devolve o valor no corte mais generoso, com a variação como incerteza.
    Se a variação não estiver encolhendo, dispara alerta: a integral pode
    simplesmente não convergir (foi o que aconteceu com rho_d em d < 3).
    """
    vals = [float(func(c)) for c in cortes]
    difs = [abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1)]
    alertas = []
    if len(difs) >= 2 and difs[-1] > difs[-2]:
        alertas.append("variação com o corte está CRESCENDO — pode não convergir")
    inc = difs[-1] if difs else 0.0
    return Medida(valor=vals[-1], incerteza=inc, unidade=unidade,
                  orcamento={"truncamento": inc}, alertas=alertas)


def combinar(*medidas, unidade=""):
    """Soma incertezas em quadratura (fontes independentes)."""
    val = medidas[-1].valor
    orc, al = {}, []
    for m in medidas:
        orc.update(m.orcamento)
        al.extend(m.alertas)
    inc = math.sqrt(sum(v**2 for k, v in orc.items() if "ordem" not in k))
    return Medida(val, inc, unidade or medidas[-1].unidade, orc, al)


# =============================================================================
#  3. Portões de validade — antes de comparar com teoria universal
# =============================================================================
LIMIAR_UNIVERSAL = 10.0     # |a|/r0 mínimo para "regime universal"


def portao_universalidade(a, r0, limiar=LIMIAR_UNIVERSAL):
    """A teoria de alcance zero vale aqui?  Precisa de |a| >> r0."""
    if r0 == 0:
        return []
    razao = abs(a) / abs(r0)
    if math.isinf(razao):
        return []
    if razao < limiar:
        return [f"|a|/r0 = {razao:.2f} < {limiar:g}: FORA do regime universal — "
                f"não compare com teoria de alcance zero"]
    return []


def portao_onda_s(k, R):
    """Só onda-s? Precisa de k*R << 1 (a onda-p entra como (kR)^4)."""
    kR = abs(k * R)
    if kR > 0.3:
        return [f"kR = {kR:.2f}: ondas parciais l>0 podem contribuir"]
    return []


def portao_confinamento(r0, l_ho):
    """Alcance zero numa armadilha exige r0 << l_ho."""
    x = abs(r0 / l_ho)
    if x > 0.1:
        return [f"r0/l_ho = {x:.3f} > 0.1: correção de alcance finito não é pequena"]
    return []


def checar(a=None, r0=None, k=None, R=None, l_ho=None):
    """Roda todos os portões que fizerem sentido com o que você passou."""
    al = []
    if a is not None and r0 is not None:
        al += portao_universalidade(a, r0)
    if k is not None and R is not None:
        al += portao_onda_s(k, R)
    if r0 is not None and l_ho is not None:
        al += portao_confinamento(r0, l_ho)
    return al


# =============================================================================
#  4. Comparação com a literatura, já com tudo junto
# =============================================================================
def inc_de_arredondamento(x_str) -> float:
    """Incerteza implicita de um valor publicado: 5.40 -> +-0.005."""
    t = str(x_str)
    return 0.5 * 10 ** (-len(t.split(".")[1])) if "." in t else 0.5


def comparar(medida: Medida, referencia: float, rotulo: str = "",
             tol_sigma: float = 3.0, inc_ref=None):
    """Linha de conformidade. inc_ref = incerteza da LITERATURA.

    Se o artigo publica 5.40, a incerteza dele e +-0.005 e ela DOMINA a nossa.
    Comparar um numero de 11 digitos com um de 3 e sempre 'FORA' — e o erro
    esta na comparacao, nao no calculo.
    """
    if inc_ref is None:
        inc_ref = inc_de_arredondamento(referencia)
    sigma_tot = math.hypot(medida.incerteza, inc_ref)
    ok = abs(medida.valor - referencia) <= tol_sigma * sigma_tot
    if sigma_tot > 0:
        n_sigma = abs(medida.valor - referencia) / sigma_tot
        sig = f"{n_sigma:5.1f}σ"
    else:
        sig = "  —  "
    status = "OK " if ok else "FORA"
    if medida.alertas:
        status = "?? "
    print(f"  {rotulo:<28} {str(medida):<34} ref={referencia:<12.6g} {sig}  {status}")
    for a in medida.alertas:
        print(f"      [!] {a}")
    return ok

# -*- coding: utf-8 -*-
"""
Benchmarks construídos a partir das REFERÊNCIAS do artigo (ver
referencias/triagem_referencias.md):

  [19] Viterbo, Lemes & Braga, RBEF 36, 1310 (2014) — equação de fase
       variável.  Em energia zero (Calogero), o "comprimento de
       espalhamento acumulado" a(r) obedece

           a'(r) = 2 V̄(r) [r - a(r)]²,   a(0) = 0,   a = a(R).

       É um método INDEPENDENTE do nosso (ODE de 1a ordem não linear,
       RK4) e serve de validação cruzada.  a(r) tem polos quando o
       potencial suporta estados ligados; integramos alternando entre
       a(r) e w(r) = 1/a(r), que obedece  w' = -2 V̄ (r w - 1)².

  [17] Pera & Boronat, Am. J. Phys. 91, 90 (2023) — a e r0 analíticos
       para esfera dura/mole (entre outros).  Esfera mole (barreira
       repulsiva V = +v μ², r < R):  a = R - tanh(κR)/κ, κ = sqrt(2v)/R.
       Limite de esfera dura (v → ∞): a → R e r0 → 2R/3.

  [22] Cencek et al., JCP 136, 224303 (2012) e [25] Janzen & Aziz,
       JCP 103 (1995) — dímero de 4He.  Aqui rodamos o modelo
       Lennard-Jones clássico do hélio (parâmetros de de Boer:
       eps/kB = 10.22 K, sigma = 2.556 Å) em unidades Å/K e comparamos
       com os valores ab initio (a = 90.4 Å, E = -1.62 mK) — a diferença
       quantifica o quanto o LJ é (ou não) um bom modelo para o He2.
"""
import csv
import math
import os

import numpy as np

from ..dois_corpos import espalhamento, analitico
from ..dois_corpos.potenciais import Potencial, PocoEsferico, PoschlTeller, Gaussiano, \
    LennardJones
from ..comum import constantes as cte

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_REF = os.path.join(AQUI, "referencias")


# ================================================== [19] fase variável
def a_fase_variavel(pot, dr=1e-4):
    """Integra a equação de fase variável de energia zero por RK4,
    alternando entre a(r) e w = 1/a perto dos polos."""
    U = lambda r: 2.0 * np.asarray(pot.V(r))          # u'' = U u
    fa = lambda r, a: float(U(r)) * (r - a)**2        # a' = U (r-a)^2
    fw = lambda r, w: -float(U(r)) * (r * w - 1.0)**2  # w' = -U (rw-1)^2

    def rk4(f, r, y, h):
        k1 = f(r, y)
        k2 = f(r + h / 2, y + h * k1 / 2)
        k3 = f(r + h / 2, y + h * k2 / 2)
        k4 = f(r + h, y + h * k3)
        return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    r = pot.r_min
    modo, y = "a", 0.0                     # começa em a(r_min) = 0
    while r < pot.R:
        h = min(dr, pot.R - r)
        y = rk4(fa if modo == "a" else fw, r, y, h)
        r += h
        # troca de representação para evitar os polos de a(r)
        if modo == "a" and abs(y) > 10.0:
            modo, y = "w", 1.0 / y
        elif modo == "w" and abs(y) > 10.0:
            modo, y = "a", 1.0 / y
    return y if modo == "a" else 1.0 / y


def bench_fase_variavel():
    """Compara a fase variável com o nosso Numerov + Eq. (110)."""
    casos = [
        ("poco nn", PocoEsferico(1.1096, 0.3918)),
        ("mpt nn", PoschlTeller(0.9071, 0.7991)),
        ("gauss nn", Gaussiano(1.2121, 0.5672)),
        ("gauss deuteron", Gaussiano(1.9102, 0.6754)),   # tem polo (1 nó)
        ("mpt deuteron", PoschlTeller(1.4388, 0.8631)),
    ]
    linhas = []
    for nome, pot in casos:
        a_fv = a_fase_variavel(pot, dr=1e-4)
        a_num = espalhamento.calcular(pot, dr=1e-3).a
        linhas.append({"benchmark": "fase_variavel[19]", "caso": nome,
                       "a_fase_var": a_fv, "a_numerov": a_num,
                       "erro_rel": abs(a_fv / a_num - 1.0)})
    return linhas


# ============================================ [17] esfera mole / dura
class EsferaMole(Potencial):
    """Barreira repulsiva: V = +v mu^2 para r <= R = 1/mu; 0 fora."""
    nome = "esfera_mole"
    rotulo = "Esfera mole"

    def __init__(self, v, mu):
        self.v = float(v)
        self.mu = float(mu)
        self.R = 1.0 / self.mu

    def V(self, r):
        r = np.asarray(r, dtype=float)
        return np.where(r <= self.R * (1.0 + 1e-12), self.v * self.mu**2, 0.0)

    def parametros(self):
        return {"v": self.v, "mu": self.mu}


def a_esfera_mole(v, R=1.0):
    """Analítico (Pera & Boronat): a = R - tanh(kR)/k, kR = sqrt(2v)."""
    k = math.sqrt(2.0 * v) / R
    return R - math.tanh(k * R) / k


def bench_esfera():
    linhas = []
    for v in (0.5, 2.0, 10.0, 1e4, 1e8):     # v -> inf: esfera dura
        pot = EsferaMole(v, 1.0)
        res = espalhamento.calcular(pot, dr=5e-4)
        a_ana = a_esfera_mole(v, 1.0)
        linhas.append({"benchmark": "esfera_mole[17]", "caso": f"v={v:g}",
                       "a_numerico": res.a, "a_analitico": a_ana,
                       "erro_rel_a": abs(res.a / a_ana - 1.0),
                       "r0_numerico": res.r0,
                       "r0_esfera_dura_2R/3": 2.0 / 3.0})
    return linhas


# ======================================= [22]/[25] dímero de 4He (LJ)
def bench_he_dimer():
    """LJ do hélio (de Boer): eps/kB = 10.22 K, sigma = 2.556 Å, em
    unidades l = 1 Å, energia eps_u = hbar²/(m_r Å²).  Compara com os
    valores ab initio de [22]: a = 90.4 Å, r0 = 8.0 Å, E = -1.62 mK."""
    eps_K, sigma = 10.22, 2.556                     # K, Å
    h2_2mr = cte.hbar2_over_2mr_K_A2(cte.M_R_HE_U)  # K.Å² (por kB)
    eps_u = 2.0 * h2_2mr                            # hbar²/(m_r Å²) em K
    # V = 4 eps [(s/r)^12 - (s/r)^6] = (1/2)(C12/r^12 - C6/r^6) * eps_u
    C12 = 8.0 * (eps_K / eps_u) * sigma**12
    C6 = 8.0 * (eps_K / eps_u) * sigma**6
    pot = LennardJones(C12, C6)
    res = espalhamento.calcular(pot, dr=2e-3)
    E_fr = analitico.energia_fr(res.a, res.r0, h2_2mr) * 1e3   # mK
    E_zr = analitico.energia_zr(res.a, h2_2mr) * 1e3           # mK
    return [{"benchmark": "he_dimer_LJ[22,25]",
             "C12_A10": C12, "C6_A4": C6,
             "r_min_A": pot.r_min, "R_A": pot.R,
             "a_A": res.a, "r0_A": res.r0, "nos": res.nos,
             "E_zr_mK": E_zr, "E_fr_mK": E_fr,
             "a_ab_initio_A": 90.4, "r0_ab_initio_A": 8.0,
             "E_ab_initio_mK": -1.62}]


# ------------------------------------------------------------- escrita
def rodar_tudo():
    os.makedirs(DIR_REF, exist_ok=True)
    blocos = [("fase_variavel", bench_fase_variavel()),
              ("esfera_mole", bench_esfera()),
              ("he_dimer", bench_he_dimer()),
              ("gauss_jeszenszki", bench_gauss_jeszenszki()),
              ("aziz", bench_aziz()),
              ("singleto_np", bench_singleto_np())]
    for nome, linhas in blocos:
        caminho = os.path.join(DIR_REF, f"bench_{nome}.csv")
        chaves = []
        for ln in linhas:
            for k in ln:
                if k not in chaves:
                    chaves.append(k)
        with open(caminho, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=chaves)
            w.writeheader()
            w.writerows(linhas)
        print("  ->", caminho)
        for ln in linhas:
            print("     ", {k: (f"{v:.6g}" if isinstance(v, float) else v)
                            for k, v in ln.items()})
    return blocos


if __name__ == "__main__":
    rodar_tudo()


# ===================================== [34] gaussiano (Jeszenszki et al.)
# PRA 97, 042708 (2018); arXiv:1802.07063.  Eles escrevem o potencial como
# V(r) = -V0 exp(-r²/L²) e dão a fórmula aproximada (Eq. 25)
#
#     a/L = soma_i alpha_i * V0 / (V0 - W_i),   i = 1..4,
#
# com W_i (posições dos polos = limiares de estados ligados) e alpha_i
# ajustados numericamente (Tabela I, linha n=4; erro relativo < 1e-4 na
# região ajustada V0 em [0, 14] (h²/µL²)).
#
# Mapeamento para a nossa convenção (u'' = 2 V̄ u, V̄ = -v µ² e^{-r²µ²}):
#     L = 1/µ  e  V0[h²/µL²] = 2 v   (o fator 2 vem do u'' = 2V̄u).
# Checagem independente: W1/2 = 1.34200232... = nosso v crítico ajustado
# na unitariedade (1.3420) — concordância externa do limiar.
W_JESZ = (2.68400465092, 17.7956995472, 45.5734799205, 85.9634003809)
ALPHA_JESZ = (1.12034897387, 0.326461774698, 0.135560767226, 0.375312300726)


def a_gauss_jeszenszki(v, mu=1.0):
    """a do potencial gaussiano pela Eq. (25) de [34] (n = 4)."""
    V0 = 2.0 * v
    soma = sum(al * V0 / (V0 - W) for al, W in zip(ALPHA_JESZ, W_JESZ))
    return soma / mu                              # a = (a/L) * L, L = 1/mu


def bench_gauss_jeszenszki():
    """Varre v e compara a numérica com a Eq. (25) de [34]."""
    linhas = []
    for v in (0.25, 0.5, 0.75, 1.0, 1.2121, 1.5, 2.0, 3.0, 5.0, 7.0):
        if abs(2 * v - W_JESZ[0]) < 0.05:          # evita o polo
            continue
        res = espalhamento.calcular(Gaussiano(v, 1.0), dr=1e-3)
        a_ref = a_gauss_jeszenszki(v, 1.0)
        linhas.append({"benchmark": "gauss_jeszenszki[34]", "v": v,
                       "a_numerico": res.a, "a_eq25": a_ref,
                       "erro_rel": abs(res.a / a_ref - 1.0),
                       "nos": res.nos})
    # limiar do 1º estado ligado: W1/2 vs nosso ajuste da unitariedade
    linhas.append({"benchmark": "gauss_limiar[34]", "v": W_JESZ[0] / 2,
                   "a_numerico": None, "a_eq25": None,
                   "erro_rel": abs(W_JESZ[0] / 2 / 1.3420 - 1.0),
                   "nos": "v_critico vs ajuste 1.3420"})
    return linhas


# ============================== [25] potencial de Aziz HFD-B do hélio
class AzizHFDB(Potencial):
    """He-He HFD-B(HE) (Aziz et al. 1987; usado em Janzen & Aziz [25]):

        V(r) = eps * [ A e^{-alpha x + beta x^2}
                       - F(x) (c6/x^6 + c8/x^8 + c10/x^10) ],  x = r/rm,
        F(x) = exp[-(D/x - 1)^2] se x < D; 1 se x >= D.

    Parâmetros (eps em K, rm em Angstrom):
        eps = 10.948, rm = 2.963, A = 1.8443101e5, alpha = 10.43329537,
        beta = -2.27965105, c6 = 1.36745214, c8 = 0.42123807,
        c10 = 0.17473318, D = 1.4826.

    Diferente do LJ, é FINITO na origem (V(0) = eps*A ~ 2e6 K), então
    r_min = 0.  Unidades do laboratório: l = 1 Angstrom,
    eps_u = hbar^2/(m_r Angstrom^2) com m_r = m(4He)/2.
    """
    nome = "aziz"
    rotulo = "Aziz HFD-B (He)"
    EPS_K = 10.948
    RM = 2.963
    A = 1.8443101e5
    ALPHA = 10.43329537
    BETA = -2.27965105
    C6, C8, C10 = 1.36745214, 0.42123807, 0.17473318
    D = 1.4826

    def __init__(self):
        self.eps_u = 2.0 * cte.hbar2_over_2mr_K_A2(cte.M_R_HE_U)  # K
        self.eps_bar = self.EPS_K / self.eps_u       # profundidade adim.
        # alcance numérico: cauda dominada por -eps c6 (rm/r)^6
        self.R = self.RM * (self.eps_bar * self.C6 / 1e-15)**(1.0 / 6)

    def _V_K(self, r):
        """V em Kelvin (forma HFD-B)."""
        x = np.asarray(r, dtype=float) / self.RM
        x = np.where(x > 0, x, 1e-12)
        F = np.where(x < self.D, np.exp(-(self.D / x - 1.0)**2), 1.0)
        rep = self.A * np.exp(-self.ALPHA * x + self.BETA * x**2)
        atr = F * (self.C6 / x**6 + self.C8 / x**8 + self.C10 / x**10)
        return self.EPS_K * (rep - atr)

    def V(self, r):
        return self._V_K(r) / self.eps_u             # adimensional

    def parametros(self):
        return {"eps_K": self.EPS_K, "rm_A": self.RM}


def bench_aziz(dr=2e-3):
    """Fecha a história do dímero de He: o LJ clássico NÃO liga
    (a = -178 A), o Aziz HFD-B liga?  Compare com ab initio [22]:
    a = 90.4 A, r0 = 8.0 A, E = -1.62 mK."""
    pot = AzizHFDB()
    res = espalhamento.calcular(pot, dr=dr)
    h22m = cte.hbar2_over_2mr_K_A2(cte.M_R_HE_U)
    E_fr = analitico.energia_fr(res.a, res.r0, h22m) * 1e3 if res.a > 0 else None
    E_zr = analitico.energia_zr(res.a, h22m) * 1e3
    return [{"benchmark": "aziz_hfdb[25]", "a_A": res.a, "r0_A": res.r0,
             "nos": res.nos, "E_zr_mK": E_zr, "E_fr_mK": E_fr,
             "a_ab_initio_A": 90.4, "r0_ab_initio_A": 8.0,
             "E_ab_initio_mK": -1.62,
             "profundidade_minimo_K": float(min(pot._V_K(
                 np.linspace(2.0, 5.0, 2000))))}]


# ================================ [23] canal singleto 1S0 do n-p
def bench_singleto_np(dr=1e-3):
    """O n-p tem DOIS canais de spin: tripleto (dêuteron, a = +5.41 fm)
    e singleto 1S0 (a_s ~ -23.74 fm, r_s ~ 2.77 fm [23]).  a < 0 =>
    NÃO existe 'dêuteron singleto'; em vez de estado ligado há um
    ESTADO VIRTUAL: a Eq. (95) dá kappa < 0, e |E| ~ 66 keV vira a
    'energia do quase-estado'.  Ajustamos um mPT a (a_s, r_s) e
    verificamos 0 nós."""
    from ..dois_corpos import ajuste as aj
    a_s, r_s = -23.74, 2.77                      # fm (valores padrão [23])
    res_aj = aj.ajustar("mpt", a_s, r_s, 0.9, 0.75, dr=dr, nos_alvo=0)
    h22m = cte.hbar2_over_2mr_MeV_fm2(cte.M_R_NP_MEV)
    kappa = (1.0 - math.sqrt(1.0 - 2.0 * r_s / a_s)) / r_s   # < 0: virtual
    E_virt = -h22m * kappa**2                     # magnitude do polo virtual
    return [{"benchmark": "singleto_np[23]", "a_alvo_fm": a_s,
             "r0_alvo_fm": r_s, "v_mPT": res_aj.p1, "mu_mPT": res_aj.p2,
             "a_obtido": res_aj.a, "r0_obtido": res_aj.r0,
             "nos": int(res_aj.nos), "convergiu": bool(res_aj.convergiu),
             "kappa_fm^-1": kappa, "E_virtual_keV": E_virt * 1e3,
             "conclusao": "kappa<0: estado VIRTUAL, nao ligado"}]

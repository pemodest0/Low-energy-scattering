# -*- coding: utf-8 -*-
"""
Gera as figuras do laboratório (reprodução das Figs. 7-10 do artigo,
mais o estudo de convergência), salvas como .png em figuras/.

Todas as figuras usam as unidades adimensionais do artigo:
distâncias em unidades de r0 (Figs. 7-8) ou fm (Fig. 9), e potencial
em unidades de hbar^2/(m_r r0^2) (Figs. 7-8).
"""
import csv
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import espalhamento, analitico
from .potenciais import LennardJones, FABRICAS
from .tabelas_artigo import TABELA3, TABELA4, ROTULOS_CASO

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_RES = os.path.join(AQUI, "resultados")

CORES = {"poco": "tab:blue", "mpt": "tab:orange",
         "gauss": "tab:green", "lj": "tab:red"}
ROTULOS = {"poco": "Poço esférico", "mpt": "Pöschl-Teller mod.",
           "gauss": "Gaussiano", "lj": "Lennard-Jones"}

plt.rcParams.update({"font.size": 11, "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.3})


def _salva(fig, nome):
    os.makedirs(DIR_FIG, exist_ok=True)
    caminho = os.path.join(DIR_FIG, nome)
    fig.savefig(caminho, bbox_inches="tight")
    plt.close(fig)
    print("  ->", caminho, flush=True)


# ------------------------------------------------------------------ Fig. 7
def fig7():
    """Poço, mPT e Gaussiana na unitariedade; eixos r/r0 e V m_r r0^2/hbar^2.
    Como os parâmetros da Tabela 3 dão r0 = 1 fm, Vbar = V e rbar = r."""
    r0 = 1.0
    r = np.linspace(1e-6, 3.0, 800)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for pot in ("poco", "mpt", "gauss"):
        p = TABELA3[("unitario", pot)]
        obj = FABRICAS[pot](p["p1"], p["p2"])
        ax.plot(r / r0, obj.V(r) * r0**2, color=CORES[pot], label=ROTULOS[pot])
    ax.set_xlabel(r"$r/r_0$")
    ax.set_ylabel(r"$\bar V = m_r r_0^2 V/\hbar^2$")
    ax.set_title("Fig. 7 — potenciais atrativos na unitariedade")
    ax.legend()
    ax.set_xlim(0, 3)
    _salva(fig, "fig7_potenciais_unitariedade.png")


# ------------------------------------------------------------------ Fig. 8
def fig8():
    """Lennard-Jones na unitariedade (caroço repulsivo + cauda atrativa)."""
    r0 = 1.0
    p = TABELA4[("unitario", "lj")]
    obj = LennardJones(p["p2"], p["p1"])        # (C12, C6)
    r = np.linspace(0.05, 3.0, 2000)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(r / r0, obj.V(r) * r0**2, color=CORES["lj"], label=ROTULOS["lj"])
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_ylim(-6, 6)
    ax.set_xlim(0, 3)
    ax.set_xlabel(r"$r/r_0$")
    ax.set_ylabel(r"$\bar V = m_r r_0^2 V/\hbar^2$")
    ax.set_title("Fig. 8 — potencial de Lennard-Jones na unitariedade")
    ax.legend()
    _salva(fig, "fig8_lennard_jones.png")


# ------------------------------------------------------------------ Fig. 9
def fig9(dr=1e-3, metodo="numerov"):
    """Soluções radiais reduzidas u(r), 4 potenciais x 3 casos, com as
    soluções analíticas tracejadas (poço: Eq. 90; mPT unitário: Eq. 118)
    e a reta assintótica g0 = 1 - r/a pontilhada."""
    casos = ["nn", "unitario", "deuteron"]
    xmax = {"nn": 8.0, "unitario": 6.0, "deuteron": 8.0}
    fig, eixos = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, caso in zip(eixos, casos):
        a_med = []
        for pot in ("poco", "mpt", "gauss", "lj"):
            tab = TABELA4 if pot == "lj" else TABELA3
            p = tab[(caso, pot)]
            res = espalhamento.calcular(FABRICAS[pot](p["p1"], p["p2"]),
                                        dr=dr, metodo=metodo)
            m = res.r <= xmax[caso]
            ax.plot(res.r[m], res.u[m], color=CORES[pot],
                    lw=1.5, label=ROTULOS[pot])
            a_med.append(res.a)
        # solução analítica do poço (Eq. 90), tracejada
        p = TABELA3[(caso, "poco")]
        R = 1.0 / p["p2"]
        rr = np.linspace(0, xmax[caso], 400)
        ax.plot(rr, analitico.u_poco(rr, p["p1"], R), "k--", lw=1.0,
                label="analítico (Eq. 90)")
        # mPT na unitariedade (Eq. 118), tracejada
        if caso == "unitario":
            p = TABELA3[(caso, "mpt")]
            ax.plot(rr, analitico.u_mpt_unitario(rr, p["p2"]), "--",
                    color="gray", lw=1.0, label="mPT analítico (Eq. 118)")
        # reta assintótica g0 = 1 - r/a
        a = a_med[0]
        g0 = 1.0 - rr / a
        if not math.isinf(a) and abs(a) < 1e4:
            ax.plot(rr, g0, ":", color="k", lw=1.0, label=r"$g_0 = 1 - r/a$")
        else:
            ax.plot(rr, np.ones_like(rr), ":", color="k", lw=1.0,
                    label=r"$g_0 = 1$")
        ax.set_xlabel(r"$r$ (fm)")
        ax.set_title(ROTULOS_CASO[caso])
        ax.set_xlim(0, xmax[caso])
    eixos[0].set_ylabel(r"$u(r)$ (normalizada: $u(R) = 1 - R/a$)")
    eixos[0].legend(fontsize=8)
    fig.suptitle("Fig. 9 — soluções radiais reduzidas de energia zero", y=1.02)
    _salva(fig, "fig9_solucoes_radiais.png")


# ----------------------------------------------------------------- Fig. 10
def fig10(dr=2e-3):
    """a/r0 em função da intensidade da interação.
    (a) v para poço/mPT/Gaussiana (mu fixo nos valores da unitariedade;
        a/r0 depende só de v por invariância de escala);
    (b) C6 para o LJ (C12 fixo no valor da unitariedade).
    Tracejado: analítico (Eq. 80 para o poço; Eq. 117 para o mPT, com
    r0 numérico)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    vs = np.linspace(0.05, 2.4, 90)
    for pot in ("poco", "mpt", "gauss"):
        mu = TABELA3[("unitario", pot)]["p2"]
        num = []
        for v in vs:
            res = espalhamento.calcular(FABRICAS[pot](v, mu), dr=dr)
            num.append(res.a / res.r0)
        ax1.plot(vs, num, color=CORES[pot], lw=1.4, label=ROTULOS[pot])
    # analíticos tracejados
    ana_poco, ana_mpt = [], []
    for v in vs:
        ana_poco.append(analitico.a_poco(v, 1.0) / analitico.r0_poco(v, 1.0))
        res = espalhamento.calcular(FABRICAS["mpt"](v, 2.0), dr=dr)
        ana_mpt.append(analitico.a_mpt(v, 2.0) / res.r0)
    ax1.plot(vs, ana_poco, "--", color="k", lw=1.0, label="poço, Eq. (80)")
    ax1.plot(vs, ana_mpt, "--", color="gray", lw=1.0, label="mPT, Eq. (117)")
    ax1.set_ylim(-12, 12)
    ax1.axhline(0, color="k", lw=0.6)
    ax1.set_xlabel(r"intensidade $v$")
    ax1.set_ylabel(r"$a/r_0$")
    ax1.set_title("(a) potenciais atrativos")
    ax1.legend(fontsize=8)

    C12 = TABELA4[("unitario", "lj")]["p2"]
    c6s = np.linspace(0.05, 1.2, 60)
    num = []
    for c6 in c6s:
        res = espalhamento.calcular(LennardJones(C12, c6), dr=dr)
        num.append(res.a / res.r0)
    ax2.plot(c6s, num, color=CORES["lj"], lw=1.4, label=ROTULOS["lj"])
    ax2.set_ylim(-12, 12)
    ax2.axhline(0, color="k", lw=0.6)
    ax2.set_xlabel(r"$C_6$ (fm$^4$)")
    ax2.set_ylabel(r"$a/r_0$")
    ax2.set_title("(b) Lennard-Jones ($C_{12}$ fixo)")
    ax2.legend(fontsize=8)

    fig.suptitle("Fig. 10 — comprimento de espalhamento vs intensidade", y=1.02)
    _salva(fig, "fig10_a_vs_intensidade.png")


# ------------------------------------------------------------ convergência
def fig_convergencia():
    """Erro relativo de a e r0 vs dr (lê resultados/convergencia.csv)."""
    caminho = os.path.join(DIR_RES, "convergencia.csv")
    with open(caminho, encoding="utf-8") as f:
        linhas = list(csv.DictReader(f))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    estilos = {"central": "o-", "numerov": "s--"}
    # (a) erro em a: poço nn e mPT nn
    for pot, cor in (("poco", CORES["poco"]), ("mpt", CORES["mpt"])):
        for met in ("central", "numerov"):
            pts = [(float(x["dr"]), float(x["erro_a"])) for x in linhas
                   if x["potencial"] == pot and x["metodo"] == met
                   and x["caso"] == "nn"]
            pts.sort()
            drs = [p[0] for p in pts]
            errs = [max(p[1], 1e-16) for p in pts]
            ax1.loglog(drs, errs, estilos[met], color=cor, ms=4,
                       label=f"{ROTULOS[pot]}, {met}")
    ax1.set_xlabel(r"$\Delta r$ (fm)")
    ax1.set_ylabel(r"erro relativo de $a$")
    ax1.set_title("(a) convergência de $a$ (caso n-n)")
    ax1.legend(fontsize=8)

    # (b) erro em r0: mPT unitário (r0 = 2/mu analítico), trapézio vs Simpson
    for met in ("central", "numerov"):
        for quad, ls in (("erro_r0_trap", ":"), ("erro_r0_simpson", "-")):
            pts = [(float(x["dr"]), float(x[quad])) for x in linhas
                   if x["potencial"] == "mpt" and x["metodo"] == met
                   and x["caso"] == "unitario" and x[quad]]
            pts.sort()
            ax2.loglog([p[0] for p in pts], [max(p[1], 1e-16) for p in pts],
                       ls, marker=estilos[met][0], ms=4,
                       label=f"{met}, {'trapézio' if 'trap' in quad else 'Simpson'}")
    ax2.set_xlabel(r"$\Delta r$ (fm)")
    ax2.set_ylabel(r"erro relativo de $r_0$")
    ax2.set_title("(b) convergência de $r_0$ (mPT unitário)")
    ax2.legend(fontsize=8)

    fig.suptitle("Estudo de convergência", y=1.02)
    _salva(fig, "fig_convergencia.png")


def gerar_todas():
    print("Fig. 7 ...", flush=True); fig7()
    print("Fig. 8 ...", flush=True); fig8()
    print("Fig. 9 ...", flush=True); fig9()
    print("Fig. 10 ...", flush=True); fig10()
    print("convergência ...", flush=True); fig_convergencia()


if __name__ == "__main__":
    gerar_todas()

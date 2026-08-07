# -*- coding: utf-8 -*-
"""
Laboratório interativo de espalhamento de baixa energia (onda-s, E = 0).

Uso:  python interativo/laboratorio.py   (na raiz do projeto)

Controles:
  - rádio esquerdo:  potencial (poço, mPT, Gaussiano, LJ) e método
    (Numerov / diferença central);
  - sliders: parâmetros do potencial (v, mu ou C6, C12) e passo dr
    (em escala log10);
  - painel (a): potencial V(r);  painel (b): u(r) com a reta de
    extrapolação g0 = 1 - r/a e o intercepto r = a marcado;
  - caixa de texto: a, r0 (trapézio e Simpson) e número de nós;
  - botões "alvo: ...": rodam o ajuste da Seção 4.5 para o caso da
    Tabela 2 escolhido e movem os sliders para os parâmetros ajustados.

Dependências: numpy e matplotlib apenas.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

# permite rodar de qualquer diretório
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import espalhamento, ajuste                      # noqa: E402
from src.potenciais import FABRICAS                       # noqa: E402
from src.tabelas_artigo import TABELA2, TABELA3, TABELA4  # noqa: E402

# ------------------------------------------------------------ configuração
POTS = ["poco", "mpt", "gauss", "lj"]
ROTULOS = {"poco": "Poço esférico", "mpt": "Pöschl-Teller",
           "gauss": "Gaussiano", "lj": "Lennard-Jones"}
# faixas dos sliders: (p1_min, p1_max, p1_ini), (p2_min, p2_max, p2_ini)
FAIXAS = {
    "poco": ((0.05, 3.0, 1.1096), (0.1, 3.0, 0.3918)),
    "mpt": ((0.05, 3.0, 0.9071), (0.1, 3.0, 0.7991)),
    "gauss": ((0.05, 3.0, 1.2121), (0.1, 3.0, 0.5672)),
    "lj": ((0.05, 12.0, 6.81472), (1e-4, 4.0, 0.90485319)),  # (C6, C12)
}
NOMES = {"poco": ("v", "µ"), "mpt": ("v", "µ"),
         "gauss": ("v", "µ"), "lj": ("C6", "C12")}
DR_INI = -2.3      # log10(dr) inicial  (dr ~ 5e-3: responsivo)

estado = {"pot": "poco", "metodo": "numerov"}

fig = plt.figure(figsize=(12.5, 7))
fig.canvas.manager.set_window_title("Laboratório de espalhamento — RBEF 45, e20230079")
ax_V = fig.add_axes([0.32, 0.56, 0.42, 0.36])
ax_u = fig.add_axes([0.32, 0.10, 0.42, 0.36])
ax_txt = fig.add_axes([0.78, 0.56, 0.20, 0.36]); ax_txt.axis("off")

# ------------------------------------------------------------- widgets
ax_radio_pot = fig.add_axes([0.02, 0.66, 0.16, 0.26])
ax_radio_pot.set_title("Potencial", fontsize=10)
radio_pot = RadioButtons(ax_radio_pot, [ROTULOS[p] for p in POTS])

ax_radio_met = fig.add_axes([0.02, 0.46, 0.16, 0.14])
ax_radio_met.set_title("Método", fontsize=10)
radio_met = RadioButtons(ax_radio_met, ["Numerov", "Dif. central"])

ax_s1 = fig.add_axes([0.06, 0.34, 0.14, 0.03])
ax_s2 = fig.add_axes([0.06, 0.28, 0.14, 0.03])
ax_s3 = fig.add_axes([0.06, 0.22, 0.14, 0.03])
s1 = Slider(ax_s1, "v", *FAIXAS["poco"][0][:2], valinit=FAIXAS["poco"][0][2])
s2 = Slider(ax_s2, "µ", *FAIXAS["poco"][1][:2], valinit=FAIXAS["poco"][1][2])
s3 = Slider(ax_s3, "log10 Δr", -3.5, -1.5, valinit=DR_INI)

# botões "reproduzir alvo"
botoes = []
for i, caso in enumerate(("nn", "unitario", "deuteron")):
    axb = fig.add_axes([0.02 + 0.062 * i, 0.10, 0.058, 0.05])
    rot = {"nn": "alvo: n-n", "unitario": "alvo: unit.", "deuteron": "alvo: dêut."}[caso]
    botoes.append((Button(axb, rot), caso))

texto = ax_txt.text(0.0, 1.0, "", va="top", family="monospace", fontsize=9)


# ------------------------------------------------------------- desenho
def desenhar(_=None):
    pot_nome = estado["pot"]
    dr = 10.0 ** s3.val
    try:
        pot = FABRICAS[pot_nome](s1.val, s2.val)
        res = espalhamento.calcular(pot, dr=dr, metodo=estado["metodo"])
    except Exception as exc:
        texto.set_text(f"erro:\n{exc}")
        fig.canvas.draw_idle()
        return

    a, r0t, r0s = res.a, res.r0_trap, res.r0_simp
    # janela de plote: região física interessante
    xmax = min(max(3.0 * (r0s if 0 < r0s < 50 else 2.0), 1.5 / max(s2.val, 0.2)), 15.0)

    ax_V.clear()
    rr = np.linspace(1e-3, xmax, 600)
    ax_V.plot(rr, pot.V(rr), color="tab:blue")
    ax_V.axhline(0, color="k", lw=0.5)
    Vmin = pot.V(rr).min()
    ax_V.set_ylim(max(Vmin * 1.15, -30), max(1.0, -0.15 * Vmin))
    ax_V.set_ylabel(r"$\bar V(r)$")
    ax_V.set_title(f"{ROTULOS[pot_nome]} — método: {estado['metodo']}", fontsize=10)
    ax_V.grid(alpha=0.3)

    ax_u.clear()
    m = res.r <= xmax
    ax_u.plot(res.r[m], res.u[m], color="tab:red", lw=1.6, label=r"$u(r)$")
    g0 = 1.0 - rr / a
    ax_u.plot(rr, g0, "k--", lw=1.0, label=r"$g_0 = 1 - r/a$")
    if 0 < a < xmax:
        ax_u.plot([a], [0.0], "ko", ms=6)
        ax_u.annotate(r"$r = a$", (a, 0), textcoords="offset points",
                      xytext=(6, 8), fontsize=9)
    ax_u.axhline(0, color="k", lw=0.5)
    ax_u.set_xlabel(r"$r$ (fm)")
    ax_u.set_ylabel(r"$u(r)$")
    ax_u.legend(fontsize=8, loc="best")
    ax_u.grid(alpha=0.3)

    n1, n2 = NOMES[pot_nome]
    a_str = f"{a:.6g}" if abs(a) < 1e5 else f"{a:.3e}"
    texto.set_text(
        f"{n1} = {s1.val:.5f}\n{n2} = {s2.val:.5f}\n"
        f"Δr = {dr:.2e} fm\nN  = {len(res.r) - 2}\n\n"
        f"a  = {a_str} fm\n"
        f"r0 (trap.)   = {r0t:.5f} fm\n"
        f"r0 (Simpson) = {r0s:.5f} fm\n\n"
        f"nós = {res.nos}\n"
        f"alcance R = {res.R_match:.3f} fm"
    )
    fig.canvas.draw_idle()


def muda_pot(rotulo):
    pot = [p for p in POTS if ROTULOS[p] == rotulo][0]
    estado["pot"] = pot
    (lo1, hi1, v1), (lo2, hi2, v2) = FAIXAS[pot]
    n1, n2 = NOMES[pot]
    s1.label.set_text(n1); s1.valmin, s1.valmax = lo1, hi1
    s1.ax.set_xlim(lo1, hi1); s1.set_val(v1)
    s2.label.set_text(n2); s2.valmin, s2.valmax = lo2, hi2
    s2.ax.set_xlim(lo2, hi2); s2.set_val(v2)
    # set_val já chama desenhar via on_changed


def muda_met(rotulo):
    estado["metodo"] = "numerov" if rotulo == "Numerov" else "central"
    desenhar()


def reproduz_alvo(caso):
    """Roda o ajuste (Seção 4.5) para o caso da Tabela 2 e move os sliders."""
    alvo = TABELA2[caso]
    pot = estado["pot"]
    tab = TABELA4 if pot == "lj" else TABELA3
    ref = tab[(caso, pot)]
    texto.set_text("ajustando...\n(aguarde)")
    fig.canvas.draw()
    try:
        aj = ajuste.ajustar(pot, alvo["a"], alvo["r0"],
                            1.1 * ref["p1"], 0.9 * ref["p2"],
                            dr=max(10.0 ** s3.val, 2e-3),
                            nos_alvo=alvo["nos"])
        s1.set_val(aj.p1)
        s2.set_val(aj.p2)
    except Exception as exc:
        texto.set_text(f"ajuste falhou:\n{exc}")
        fig.canvas.draw_idle()


radio_pot.on_clicked(muda_pot)
radio_met.on_clicked(muda_met)
s1.on_changed(desenhar)
s2.on_changed(desenhar)
s3.on_changed(desenhar)
for botao, caso in botoes:
    botao.on_clicked(lambda _, c=caso: reproduz_alvo(c))

desenhar()

if __name__ == "__main__":
    plt.show()

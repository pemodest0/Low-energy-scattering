# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

import comum
from src import espalhamento, solvers
from src.potenciais import FABRICAS

comum.setup("Estação 1 — Espalhamento de baixa energia", "🌊")
sim, teoria, codigo, resultados = comum.abas()

POTS = {"Poço esférico": "poco", "Pöschl-Teller": "mpt",
        "Gaussiano": "gauss", "Lennard-Jones": "lj"}
PRESETS = {
    "n–n (a<0)": {"poco": (1.1096, 0.3918), "mpt": (0.9071, 0.7991),
                  "gauss": (1.2121, 0.5672), "lj": (9.86668911, 3.08836698)},
    "unitariedade": {"poco": (1.2337, 1.0), "mpt": (1.0, 2.0),
                     "gauss": (1.3420, 1.4349), "lj": (0.26462461, 0.00034068)},
    "dêuteron (a>0)": {"poco": (1.7575, 0.5), "mpt": (1.4388, 0.8631),
                       "gauss": (1.9102, 0.6754), "lj": (6.81472, 0.90485319)},
}

with sim:
    c1, c2, c3 = st.columns([1.2, 1, 1])
    pot_rotulo = c1.selectbox("potencial", list(POTS))
    metodo = c2.selectbox("método", ["numerov", "central"])
    preset = c3.selectbox("preset", ["livre"] + list(PRESETS))
    nome = POTS[pot_rotulo]
    n1, n2 = ("C₆", "C₁₂") if nome == "lj" else ("v (profundidade)", "µ (1/alcance)")
    if preset != "livre":
        p1_ini, p2_ini = PRESETS[preset][nome]
    else:
        p1_ini, p2_ini = (1.0, 1.0) if nome != "lj" else (6.8, 0.9)
    c1, c2, c3 = st.columns(3)
    p1 = c1.number_input(n1, 1e-5, 50.0, float(p1_ini), format="%.5f")
    p2 = c2.number_input(n2, 1e-5, 50.0, float(p2_ini), format="%.5f")
    dr = c3.select_slider("passo Δr (fm)", [4e-3, 2e-3, 1e-3, 5e-4], 1e-3)

    pot = FABRICAS[nome](p1, p2)
    n_passos = int((pot.R - pot.r_min) / dr)
    st.caption(f"⏱️ tempo estimado: {comum.tempo_estimado(n_passos)} "
               f"(alcance R = {pot.R:.1f} fm → {n_passos:,} passos)".replace(",", "."))
    if st.button("▶️ Rodar", type="primary"):
        res = comum.rodar_medido(
            lambda: espalhamento.calcular(pot, dr=dr, metodo=metodo), n_passos)
        m1, m2, m3, m4 = st.columns(4)
        a_txt = f"{res.a:.4f}" if abs(res.a) < 1e5 else f"{res.a:.2e}"
        m1.metric("a (fm)", a_txt)
        m2.metric("r₀ Simpson (fm)", f"{res.r0_simp:.5f}")
        m3.metric("r₀ trapézio (fm)", f"{res.r0_trap:.5f}")
        m4.metric("nós", res.nos)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        xmax = min(3 * max(res.r0_simp, 1.0) + 2, res.R_match)
        rr = np.linspace(max(pot.r_min, 1e-3), xmax, 500)
        ax1.plot(rr, pot.V(rr), "b-"); ax1.axhline(0, color="k", lw=0.5)
        Vmin = float(pot.V(rr).min())
        ax1.set_ylim(max(Vmin * 1.2, -30), max(1, -0.2 * Vmin))
        ax1.set_ylabel("V̄(r)"); ax1.grid(alpha=0.3)
        m = res.r <= xmax
        ax2.plot(res.r[m], res.u[m], "r-", lw=2, label="u(r)")
        g0 = 1 - rr / res.a
        ax2.plot(rr, g0, "k--", label="g₀ = 1 − r/a")
        if 0 < res.a < xmax:
            ax2.plot([res.a], [0], "ko")
        ax2.axhline(0, color="k", lw=0.5); ax2.grid(alpha=0.3)
        ax2.set_xlabel("r (fm)"); ax2.set_ylabel("u(r)"); ax2.legend()
        st.pyplot(fig)

with teoria:
    st.markdown("**Uma frase:** fora do alcance do potencial, a solução de "
                "energia zero é a reta $g_0 = 1 - r/a$; o intercepto é o "
                "comprimento de espalhamento.")
    st.latex(r"u''(r) = 2\,\bar V(r)\,u(r), \qquad u(0)=0")
    st.latex(r"a = R - \frac{2\Delta r\, u(R)}{u(R+\Delta r)-u(R-\Delta r)}"
             r"\quad\text{(Eq. 110)}")
    st.latex(r"r_0 = 2\int_0^R\left[g_0^2 - u_0^2\right]dr \quad\text{(Eq. 56)}")
    st.markdown("Sinais: $a<0$ quase liga (n–n) · $|a|\\to\\infty$ limiar "
                "(unitariedade, regime de Efimov) · $a>0$ estado ligado raso "
                "(dêuteron, 1 nó). **Fundo:** notebooks 01–02 (`estudar.bat`).")

with codigo:
    comum.mostrar_fonte(espalhamento.calcular, solvers.integrar)

with resultados:
    comum.mostra_csv("resultados/reproducao_tabelas_3_4.csv",
                     "Parâmetros publicados no artigo → (a, r₀) do nosso código")
    comum.mostra_csv("resultados/validacao_analitica.csv",
                     "Validação contra fórmulas fechadas (poço e mPT)")

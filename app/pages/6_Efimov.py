# -*- coding: utf-8 -*-
import math

import matplotlib.pyplot as plt
import streamlit as st

import comum
from src import efimov

comum.setup("Estação 6 — Efimov: a torre de trímeros", "🗼")
sim, teoria, codigo, resultados = comum.abas()

with sim:
    st.metric("s₀ resolvido do zero (eq. transcendental)",
              f"{efimov.S0:.10f}", help="literatura: 1.00624")
    n_niveis = st.slider("níveis da torre", 2, 6, 4)
    st.caption(f"⏱️ {comum.tempo_estimado(9000 * 70 * n_niveis)}")
    if st.button("▶️ Calcular a torre", type="primary"):
        ks = comum.rodar_medido(lambda: efimov.torre(n_niveis),
                                9000 * 70 * n_niveis)
        rz = efimov.razoes(ks)
        fig, ax = plt.subplots(figsize=(7.5, 4.6))
        for n, k in enumerate(ks):
            ax.hlines(-(k**2), 0, 1 / k, lw=3, color="tab:purple")
            ax.annotate(f"  n={n}", (1 / k, -(k**2)), fontsize=9)
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-15)
        ax.set_xlabel("extensão ~ 1/κ (unid. R₀, log)")
        ax.set_ylabel("E = −κ² (log)"); ax.grid(alpha=0.3)
        st.pyplot(fig)
        alvo = math.exp(2 * math.pi / efimov.S0)
        comum.cartao_validacao(
            [(f"E{n}/E{n+1}", f"{r*r:.2f}", f"{alvo:.2f}",
              "e^(2π/s₀) universal") for n, r in enumerate(rz)][:4])
        st.caption("O 1º degrau desvia um pouco (efeito do corte de curto "
                   "alcance) — mesmo desvio que o experimento vê no ³⁹K.")

with teoria:
    st.markdown("**Uma frase:** com |a|=∞ em 2 corpos, 3 bósons ganham uma "
                "torre INFINITA de estados ligados em escala geométrica — "
                "cada nível 515× mais fraco e 22,7× maior.")
    st.latex(r"g''(x) + \left[s_0^2 - (\kappa R_0)^2 e^{2x}\right] g = 0,"
             r"\qquad x = \ln(R/R_0)")
    st.latex(r"\frac{E_{n+1}}{E_n} = e^{-2\pi/s_0} \approx \frac{1}{515{,}03}")
    st.markdown("A dissertação pergunta: **como essa torre morre quando o "
                "gás é achatado 3D→2D** (estação 4)? Em 2D puro não há "
                "Efimov. **Fundo:** notebook 03; Braaten & Hammer (2006).")

with codigo:
    comum.mostrar_fonte(efimov.s0_universal, efimov.torre)

with resultados:
    st.markdown("| checagem | nosso valor | universal |\n|---|---|---|\n"
                "| s₀ | 1,0062378 | 1,00624 |\n"
                "| razão de energias (níveis rasos) | 515,04 | 515,03 |\n"
                "| razão de tamanhos | 22,694 | e^(π/s₀) = 22,694 |")

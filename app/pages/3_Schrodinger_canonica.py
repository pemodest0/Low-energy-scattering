# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

import comum
from src import schrodinger as sq

comum.setup("Estação 3 — Schrödinger canônica", "📐")
sim, teoria, codigo, resultados = comum.abas()

with sim:
    sistema = st.radio("sistema", ["Oscilador harmônico", "Hidrogênio (radial)",
                                   "Morse (molécula)", "Tunelamento (barreira)"],
                       horizontal=True)
    if sistema == "Oscilador harmônico":
        n_est = st.slider("nº de estados", 2, 8, 4)
        st.caption(f"⏱️ {comum.tempo_estimado(4000 * 45 * n_est)}")
        if st.button("▶️ Rodar", type="primary"):
            niveis = comum.rodar_medido(
                lambda: sq.autovalores(sq.V_oscilador, -9, 9, 0.0,
                                       n_est + 1.0, n_estados=n_est, dx=4e-3),
                4000 * 45 * n_est)
            fig, ax = plt.subplots(figsize=(7.5, 5))
            xx = np.linspace(-5.5, 5.5, 300)
            ax.plot(xx, sq.V_oscilador(xx), "k-", lw=2)
            for E, n in niveis:
                x, u = sq.autofuncao(sq.V_oscilador, E, -9, 9, dx=4e-3)
                m = np.abs(x) <= 5.5
                ax.plot(x[m], E + 0.9 * u[m], lw=1.4)
                ax.axhline(E, color="gray", lw=0.5, ls=":")
            ax.set_ylim(0, n_est + 1); ax.grid(alpha=0.3)
            ax.set_xlabel("x"); ax.set_ylabel("E  /  E + ψ(x)")
            st.pyplot(fig)
            comum.cartao_validacao(
                [(f"E{n}", f"{E:.6f}", f"{n + 0.5}", "n+1/2 exato")
                 for E, n in niveis[:4]])
    elif sistema == "Hidrogênio (radial)":
        st.caption(f"⏱️ {comum.tempo_estimado(30000 * 45 * 3)}")
        if st.button("▶️ Rodar", type="primary"):
            niveis = comum.rodar_medido(
                lambda: sq.autovalores(sq.V_hidrogenio(0), 1e-6, 60,
                                       -0.6, -0.03, n_estados=3, dx=2e-3),
                30000 * 45 * 3)
            comum.cartao_validacao(
                [(f"n={nos+1}", f"{E:.6f}", f"{-0.5/(nos+1)**2:.6f}",
                  "-1/2n² exato") for E, nos in niveis])
            st.caption("Unidades atômicas: 1 u.a. = 27,2 eV. Os níveis se "
                       "acumulam em E=0 — acima disso é o CONTÍNUO, o mundo "
                       "da estação 1.")
    elif sistema == "Morse (molécula)":
        D = st.slider("profundidade D", 4.0, 20.0, 10.0)
        st.caption(f"⏱️ {comum.tempo_estimado(3000 * 45 * 4)} · "
                   f"níveis que cabem: {sq.n_max_morse(D, 1.0) + 1}")
        if st.button("▶️ Rodar", type="primary"):
            nmax = min(sq.n_max_morse(D, 1.0) + 1, 6)
            niveis = comum.rodar_medido(
                lambda: sq.autovalores(sq.V_morse(D, 1.0, 2.0), 0.05, 12,
                                       -D + 0.01, -0.02, n_estados=nmax,
                                       dx=4e-3), 3000 * 45 * nmax)
            comum.cartao_validacao(
                [(f"E{n}", f"{E:.4f}", f"{sq.E_morse(n, D, 1.0):.4f}",
                  "fórmula de Morse") for E, n in niveis[:4]])
            st.caption("Os níveis se APERTAM subindo (anarmonicidade) e o "
                       "poço comporta um número FINITO — como o He₂ real.")
    else:
        V0 = st.slider("altura da barreira V₀", 2.0, 10.0, 5.0)
        st.caption(f"⏱️ {comum.tempo_estimado(20000 * 60)}")
        if st.button("▶️ Rodar", type="primary"):
            Es = np.linspace(0.2, 2.8 * V0, 60)
            Ts = comum.rodar_medido(
                lambda: [sq.transmissao(sq.V_barreira(V0, 1.0), E, -0.5, 1.5,
                                        dx=1e-4) for E in Es], 20000 * 60)
            Ta = [sq.T_barreira_analitico(E, V0, 1.0) for E in Es]
            fig, ax = plt.subplots(figsize=(7.5, 4.2))
            ax.semilogy(Es, Ts, "r-", lw=2, label="numérico")
            ax.semilogy(Es, Ta, "k--", label="exato")
            ax.axvline(V0, color="gray", ls=":")
            ax.set_xlabel("E"); ax.set_ylabel("transmissão T (log)")
            ax.legend(); ax.grid(alpha=0.3)
            st.pyplot(fig)
            st.caption("Abaixo de V₀: queda exponencial (tunelamento). Acima: "
                       "T < 1 com ressonâncias de transmissão.")

with teoria:
    st.markdown("**Uma frase:** o mesmo Numerov do espalhamento, agora "
                "caçando as energias em que $u$ morre nas duas pontas.")
    st.latex(r"u''(x) = 2\,[V(x) - E]\,u(x), \qquad u(x_0)=u(x_1)=0")
    st.markdown("O nº de nós conta os autovalores abaixo de E (teorema da "
                "oscilação) → bisseção nos nós + `brentq`. **Fundo:** "
                "notebook 04.")

with codigo:
    comum.mostrar_fonte(sq.autovalor_n, sq.transmissao)

with resultados:
    st.markdown("| sistema | precisão vs exato |\n|---|---|\n"
                "| oscilador | 10⁻⁸ |\n| hidrogênio | 10⁻⁶ |\n"
                "| Morse | 5 casas |\n| barreira | ~10⁻⁵ |")
    st.caption("Testes automáticos correspondentes: tests/test_lab.py "
               "(estação 7 roda todos).")

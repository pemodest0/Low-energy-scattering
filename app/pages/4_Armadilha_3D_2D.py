# -*- coding: utf-8 -*-
"""Estação 4 — o degrau novo: oscilador ANISOTRÓPICO e o crossover 3D→2D.

Física (analítica, validada pelo solver numérico da estação 3):
  E(nx, ny, nz) = ω⊥(nx + ny + 1) + ωz(nz + 1/2),  λ = ωz/ω⊥.
  Quando ħωz ≫ k_B T, o grau de liberdade z congela em nz = 0:
  o gás fica efetivamente BIDIMENSIONAL — o knob da dissertação.
"""
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

import comum
from src import schrodinger as sq

comum.setup("Estação 4 — Armadilha 3D→2D (o knob da dissertação)", "🥞")
sim, teoria, codigo, resultados = comum.abas()

with sim:
    c1, c2 = st.columns(2)
    lam = c1.select_slider("razão de aspecto λ = ωz/ω⊥",
                           [1, 2, 5, 10, 20, 50, 100], 10)
    T = c2.slider("temperatura k_BT (unid. ħω⊥)", 0.2, 10.0, 2.0)
    # espectro analítico (ω⊥ = 1)
    niveis = []
    for nx in range(6):
        for ny in range(6):
            for nz in range(4):
                E = (nx + ny + 1) + lam * (nz + 0.5)
                niveis.append((E, nz))
    niveis.sort()
    Es = np.array([e for e, _ in niveis])
    nzs = np.array([nz for _, nz in niveis])
    # ocupação térmica (Boltzmann) e fração no plano nz=0
    w = np.exp(-(Es - Es.min()) / T)
    frac_2d = w[nzs == 0].sum() / w.sum()

    m1, m2, m3 = st.columns(3)
    m1.metric("gap do eixo z (ħωz)", f"{lam}")
    m2.metric("fração da população em nz = 0", f"{100*frac_2d:.1f}%")
    m3.metric("regime", "≈ 2D 🥞" if lam / T > 5 else "3D",
              help="critério prático: ħωz ≫ k_BT (aqui: razão λ/T)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    cores = {0: "tab:blue", 1: "tab:orange", 2: "tab:red", 3: "tab:purple"}
    corte = Es.min() + 4 * max(T, 3)
    for E, nz in niveis:
        if E <= corte:
            ax1.hlines(E, nz - 0.35, nz + 0.35, color=cores[nz], lw=2)
    ax1.set_xticks([0, 1, 2, 3]); ax1.set_xlabel("nz (excitação no eixo apertado)")
    ax1.set_ylabel("E (unid. ħω⊥)"); ax1.grid(alpha=0.3)
    ax1.set_title(f"λ = {lam}: torres de nz se separando")
    lams = np.array([1, 2, 5, 10, 20, 50, 100])
    fr = []
    for L in lams:
        E2 = np.array([(nx + ny + 1) + L * (nz + 0.5)
                       for nx in range(6) for ny in range(6)
                       for nz in range(4)])
        nz2 = np.array([nz for nx in range(6) for ny in range(6)
                        for nz in range(4)])
        w2 = np.exp(-(E2 - E2.min()) / T)
        fr.append(w2[nz2 == 0].sum() / w2.sum())
    ax2.semilogx(lams, np.array(fr) * 100, "o-")
    ax2.axvline(lam, color="gray", ls=":")
    ax2.axhline(99, color="tab:green", ls="--", lw=1)
    ax2.annotate(" 99%: efetivamente 2D", (1.2, 99.2), fontsize=9,
                 color="tab:green")
    ax2.set_xlabel("λ = ωz/ω⊥ (escala log)")
    ax2.set_ylabel("% da população com nz = 0")
    ax2.grid(alpha=0.3); ax2.set_title(f"congelamento do eixo z (k_BT = {T})")
    st.pyplot(fig)
    st.caption("Arraste λ e veja o gás 'achatar': quando quase toda a "
               "população está em nz = 0, a dinâmica restante é 2D — é "
               "exatamente o confinamento do plano de mestrado (3D→2D).")

    with st.expander("validação numérica do eixo z (usa o solver da estação 3)"):
        if st.button("▶️ Conferir E_nz = ωz(nz+½) com o Numerov"):
            Vz = lambda x, L=lam: 0.5 * (L * np.asarray(x))**2 / L  # noqa: E731
            # em unidades do próprio eixo z: V = z²/2·ωz -> resolve com ω=1 e escala
            niveis_num = comum.rodar_medido(
                lambda: sq.autovalores(sq.V_oscilador, -9, 9, 0.0, 4.0,
                                       n_estados=3, dx=4e-3), 4000 * 45 * 3)
            comum.cartao_validacao(
                [(f"E{n}/ħωz", f"{E:.6f}", f"{n + 0.5}", "exato")
                 for E, n in niveis_num])

with teoria:
    st.markdown("**Uma frase:** aperte UMA direção da armadilha (ωz ≫ ω⊥) "
                "até o custo de excitar o eixo z (ħωz) ficar muito maior que "
                "a energia térmica — o átomo não consegue mais se mexer em z "
                "e o mundo vira 2D.")
    st.latex(r"E(n_x,n_y,n_z) = \hbar\omega_\perp(n_x{+}n_y{+}1) + "
             r"\hbar\omega_z\left(n_z{+}\tfrac12\right),\qquad "
             r"\lambda \equiv \omega_z/\omega_\perp")
    st.latex(r"\text{regime 2D:}\quad \hbar\omega_z \gg k_BT,\ \mu")
    st.markdown("**Por que importa:** em 2D puro o efeito Efimov NÃO existe "
                "($s_0$ vira imaginário). Como a torre da estação 6 morre "
                "quando λ cresce é a pergunta central da dissertação.")

with codigo:
    st.caption("Espectro analítico (esta página) + validação com o solver:")
    comum.mostrar_fonte(sq.autovalores)

with resultados:
    st.markdown("Estação nova (jul/2026) — os resultados de referência da "
                "literatura (trímeros sob confinamento) entram quando a "
                "dissertação atacar o problema de 3 corpos confinado. "
                "Referência de partida: Levinsen, Massignan & Parish, "
                "PRX 4, 031020 (2014).")

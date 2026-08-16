# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

import comum
from src import feshbach as fb

comum.setup("Estação 5 — Feshbach do ³⁹K (o botão real)", "🧲")
sim, teoria, codigo, resultados = comum.abas()

with sim:
    estado = st.radio("estado hiperfino", list(fb.RESSONANCIAS_39K),
                      horizontal=True, format_func=lambda s: f"|1,{s.split(',')[1]}⟩")
    ressons = fb.RESSONANCIAS_39K[estado]
    idx = st.selectbox("ressonância", range(len(ressons)),
                       format_func=lambda i: f"B₀ = {ressons[i]['B0']} G "
                       f"(Δ = {ressons[i]['Delta']} G)")
    r = ressons[idx]
    span = max(abs(r["Delta"]) * 3, 30)
    B = st.slider("campo magnético B (G)",
                  float(r["B0"] - span), float(r["B0"] + span),
                  float(r["B0"] + r["Delta"] / 2))
    aB = fb.a_de_B(B, **r)
    m1, m2, m3 = st.columns(3)
    m1.metric("a(B)", f"{aB:.1f} a₀")
    m2.metric("zero de a", f"{fb.zero_de_a(r['B0'], r['Delta']):.1f} G")
    m3.metric("l_vdW do K", f"{fb.comprimento_vdw_a0():.1f} a₀")

    Bs = np.linspace(r["B0"] - span, r["B0"] + span, 1500)
    a = fb.a_de_B(Bs, **r)
    a[np.abs(Bs - r["B0"]) < span / 200] = np.nan
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.plot(Bs, a, "b-")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(r["B0"], color="tab:purple", ls="--")
    ax.plot([B], [aB], "ro", ms=9)
    ax.set_ylim(-800, 800)
    ax.set_xlabel("B (gauss)"); ax.set_ylabel("a (a₀)"); ax.grid(alpha=0.3)
    if estado == "1,1" and abs(r["B0"] - 402.5) < 1:
        Btri = fb.B_de_a(-1500.0, **r)
        ax.axvline(Btri, color="tab:red", ls=":")
        ax.annotate(" 1º trímero de Efimov\n (a₁₋ = −1500 a₀)", (Btri, -600),
                    color="tab:red", fontsize=9)
    st.pyplot(fig)
    # energia do dímero universal no lado a>0
    if aB > 50:
        lvdw = fb.comprimento_vdw_a0()
        st.caption(f"lado a>0: existe um dímero raso com E_b ≈ −ħ²/(m a²) "
                   f"(universal, válido pois a = {aB:.0f} a₀ ≫ l_vdW = "
                   f"{lvdw:.0f} a₀)")

with teoria:
    st.markdown("**Uma frase:** perto de uma ressonância de Feshbach o campo "
                "magnético vira o 'botão de v' — o mesmo a que sintonizamos "
                "nos potenciais-modelo, agora girado por uma bobina.")
    st.latex(r"a(B) = a_{bg}\left(1 - \frac{\Delta}{B - B_0}\right)")
    st.markdown("Dados: D'Errico et al., NJP 9, 223 (2007); ressonância larga "
                "do |1,1⟩ refinada em Zaccanti et al. (2009): B₀ = 402,50 G, "
                "Δ = −52,1 G. **Fundo:** notebook 03.")

with codigo:
    comum.mostrar_fonte(fb.a_de_B, fb.B_de_a, fb.comprimento_vdw_a0)

with resultados:
    st.markdown("| checagem | valor | referência |\n|---|---|---|\n"
                "| zero de a (|1,1⟩ larga) | 350,4 G | citado no artigo |\n"
                "| l_vdW do ³⁹K | 64,6 a₀ | Chin et al. 2010 |\n"
                "| B do 1º trímero | ≈ 403,5 G | a₁₋ = −1500 a₀ "
                "(Zaccanti 2009) |")

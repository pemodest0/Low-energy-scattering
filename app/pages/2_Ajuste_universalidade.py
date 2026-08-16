# -*- coding: utf-8 -*-
import math

import streamlit as st

import comum
from src import ajuste
from src.tabelas_artigo import TABELA2, TABELA3, TABELA4, NOMES_PARAM

comum.setup("Estação 2 — Ajuste e universalidade", "🎛️")
sim, teoria, codigo, resultados = comum.abas()

with sim:
    st.markdown("Sintonize um potencial para reproduzir $(a, r_0)$ alvo — "
                "o algoritmo da Seção 4.5 do artigo.")
    c1, c2 = st.columns(2)
    caso = c1.selectbox("caso físico (Tabela 2)",
                        ["nn", "unitario", "deuteron"],
                        format_func={"nn": "nêutron–nêutron",
                                     "unitario": "unitariedade",
                                     "deuteron": "dêuteron"}.get)
    pot = c2.selectbox("potencial", ["poco", "mpt", "gauss", "lj"])
    alvo = TABELA2[caso]
    ref = (TABELA4 if pot == "lj" else TABELA3)[(caso, pot)]
    n_exec = 60 if pot != "lj" else 80
    passos_por_solve = 25_000 if pot != "lj" else 200_000
    st.caption(f"⏱️ tempo estimado: "
               f"{comum.tempo_estimado(passos_por_solve, n_exec)} "
               f"(~{n_exec} soluções da equação radial)")
    if st.button("▶️ Rodar o ajuste", type="primary"):
        aj = comum.rodar_medido(
            lambda: ajuste.ajustar(pot, alvo["a"], alvo["r0"],
                                   1.1 * ref["p1"], 0.9 * ref["p2"],
                                   dr=2e-3, nos_alvo=alvo["nos"]),
            passos_por_solve * n_exec)
        n1, n2 = NOMES_PARAM[pot]
        comum.cartao_validacao([
            (f"{n1} ajustado", f"{aj.p1:.5f}", f"{ref['p1']}", "artigo"),
            (f"{n2} ajustado", f"{aj.p2:.5f}", f"{ref['p2']}", "artigo"),
            ("a obtido (fm)", f"{aj.a:.4g}",
             "±∞" if math.isinf(alvo["a"]) else alvo["a"], "Tabela 2"),
            ("r₀ obtido (fm)", f"{aj.r0:.4f}", alvo["r0"], "Tabela 2"),
        ])
        st.success(f"nós = {aj.nos} (esperado {alvo['nos']}) · "
                   f"convergiu = {aj.convergiu}")

with teoria:
    st.markdown("**Uma frase:** dois laços aninhados — a intensidade acerta "
                "$a$, o alcance acerta $r_0$; perto da unitariedade "
                "resolvemos $1/a$ (contínuo) em vez de $a$ (diverge).")
    st.latex(r"\frac{1}{a}(v) = \frac{1}{a_{\rm alvo}} \;\;\text{(laço interno)}"
             r"\qquad r_0(\mu) = r_0^{\rm alvo} \;\;\text{(laço externo)}")
    st.markdown("**Universalidade:** 4 formas de potencial completamente "
                "diferentes reproduzem os mesmos $(a, r_0)$ — a física de "
                "baixa energia não enxerga a forma. Ver Fig. 9/10 em "
                "`figuras/`.")

with codigo:
    comum.mostrar_fonte(ajuste.ajustar)

with resultados:
    comum.mostra_csv("resultados/ajustes.csv",
                     "Nossos 12 ajustes vs parâmetros publicados (~0,1–1%)")

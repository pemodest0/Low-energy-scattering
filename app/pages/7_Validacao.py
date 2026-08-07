# -*- coding: utf-8 -*-
import subprocess
import sys

import streamlit as st

import comum

comum.setup("Estação 7 — Validação (tudo vs analítico e artigo)", "✅")
sim, teoria, codigo, resultados = comum.abas()

with sim:
    st.markdown("Roda a suite completa de testes (25) — a rede de segurança "
                "do projeto.")
    st.caption("⏱️ ~20–40 s")
    if st.button("▶️ Rodar os 25 testes", type="primary"):
        with st.spinner("pytest rodando..."):
            proc = subprocess.run(
                [sys.executable, "-m", "pytest", "tests/", "-q",
                 "--no-header"], cwd=comum.RAIZ,
                capture_output=True, text=True, timeout=300)
        saida = proc.stdout.strip().splitlines()
        if proc.returncode == 0:
            st.success(saida[-1] if saida else "todos passaram")
        else:
            st.error("algum teste falhou — veja abaixo")
        st.code("\n".join(saida[-25:]), language="text")

with teoria:
    st.markdown("""**As 4 camadas de validação do projeto:**

1. **Analítico** — fórmulas fechadas (poço, mPT, esfera mole, oscilador,
   hidrogênio, Morse, barreira): erros 10⁻⁵–10⁻¹⁰.
2. **Artigo** — Tabelas 1–4 reproduzidas em ~0,1% (LJ ~1%).
3. **Métodos cruzados** — central × Numerov × fase variável: 10⁻⁶–10⁻¹³.
4. **Externo** — Jeszenszki (gaussiano), Aziz (He₂), d'Errico/Zaccanti
   (³⁹K), fator de Efimov 515,03.
""")

with codigo:
    st.caption("A suite inteira: `tests/test_lab.py` — cada teste ancora "
               "num resultado analítico ou publicado, nunca em 'valor que o "
               "código dava antes'.")

with resultados:
    for rel, tit in [
        ("resultados/validacao_analitica.csv", "Validação analítica"),
        ("resultados/reproducao_tabelas_3_4.csv", "Tabelas 3–4 do artigo"),
        ("resultados/ajustes.csv", "Ajustes vs artigo"),
        ("resultados/energias_tabela1.csv", "Energias (Tabela 1)"),
        ("referencias/bench_fase_variavel.csv", "Fase variável [19]"),
        ("referencias/bench_gauss_jeszenszki.csv", "Gaussiano [34]"),
        ("referencias/bench_aziz.csv", "Aziz He₂ [25]"),
        ("referencias/bench_he_dimer.csv", "He₂ com LJ [22]"),
    ]:
        comum.mostra_csv(rel, tit)

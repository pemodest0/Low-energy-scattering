# -*- coding: utf-8 -*-
"""Página inicial do app unificado.  Rode:  streamlit run app/Inicio.py"""
import os
import re
import sys

import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app import comum  # noqa: E402

st.set_page_config(page_title="Laboratório de Átomos Frios",
                   page_icon="🧊", layout="wide")
st.title("🧊 Laboratório de Átomos Frios — ambiente unificado")
st.caption("Teoria + código + simulação + resultados, uma estação por tema. "
           "Toda estação tem as mesmas 4 abas: Simular · Teoria · Código · Resultados.")

# métricas do projeto
try:
    s = open(os.path.join(comum.RAIZ, "pyproject.toml"), encoding="utf-8").read()
    versao = re.search(r'version\s*=\s*"([^"]+)"', s).group(1)
except OSError:
    versao = "?"
c1, c2, c3 = st.columns(3)
c1.metric("Versão", f"v{versao}")
c2.metric("Testes automáticos", "25 (estação 7 roda todos)")
c3.metric("Fonte da física", "src/ — nada duplicado aqui")

st.markdown("""
### As estações (ordem sugerida)

| # | Estação | Uma frase |
|---|---|---|
| 1 | **Espalhamento** | o laboratório original: a e r₀ de 4 potenciais, ao vivo |
| 2 | **Ajuste e universalidade** | sintonize um potencial para (a, r₀) alvo, como no artigo |
| 3 | **Schrödinger canônica** | oscilador, hidrogênio, Morse e tunelamento |
| 4 | **Armadilha 3D→2D** | o knob da dissertação: congele uma dimensão |
| 5 | **Feshbach ³⁹K** | o botão real do laboratório da Patrícia: a(B) |
| 6 | **Efimov** | a torre de trímeros com o fator 22,7 emergindo |
| 7 | **Validação** | todos os números vs analítico e artigo; rode os testes |

**Profundidade:** cada estação tem a teoria em meia página; o mergulho
completo mora nos notebooks (`estudar.bat`) e nos guias (`INICIO.html`).
""")
st.info("👈 Navegue pelas estações na barra lateral. Simulações só rodam "
        "quando você clica **Rodar** — com tempo estimado antes.")

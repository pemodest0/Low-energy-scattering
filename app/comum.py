# -*- coding: utf-8 -*-
"""
Helpers compartilhados por todas as estações do app.

Regra de ouro do app: TODA estação tem as mesmas 4 abas
(Simular | Teoria | Código | Resultados) e nenhuma física própria —
tudo vem de src/ (se o código muda, o app muda junto).
"""
import inspect
import os
import sys
import time

import streamlit as st

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)


def setup(titulo, icone="🔬"):
    st.set_page_config(page_title=titulo, page_icon=icone, layout="wide")
    st.title(f"{icone} {titulo}")


def abas():
    return st.tabs(["▶️ Simular", "📖 Teoria", "💻 Código", "📊 Resultados"])


# ------------------------------------------------ tempo estimado
def velocidade_passos():
    """Passos de Numerov por segundo (calibrado na 1ª execução)."""
    return st.session_state.get("passos_por_s", 1.5e6)


def calibra(n_passos, segundos):
    if segundos > 0.05:
        st.session_state["passos_por_s"] = n_passos / segundos


def tempo_estimado(n_passos, n_execucoes=1):
    t = n_passos * n_execucoes / velocidade_passos()
    if t < 1:
        return f"~{t:.1f} s"
    if t < 90:
        return f"~{t:.0f} s"
    return f"~{t/60:.1f} min"


def rodar_medido(fn, n_passos):
    """Executa fn() com spinner + cronômetro e calibra a velocidade."""
    t0 = time.time()
    with st.spinner("integrando..."):
        res = fn()
    dt = time.time() - t0
    calibra(n_passos, dt)
    st.caption(f"⏱️ levou {dt:.2f} s ({n_passos:,} passos)".replace(",", "."))
    return res


# ------------------------------------------------ cartões padrão
def cartao_validacao(linhas):
    """linhas: [(nome, valor_num, valor_ref, fonte)] -> métricas + erro."""
    cols = st.columns(len(linhas))
    for c, (nome, num, ref, fonte) in zip(cols, linhas):
        with c:
            st.metric(nome, f"{num}", help=f"referência: {ref} ({fonte})")
            try:
                err = abs(float(num) / float(ref) - 1.0)
                cor = "✅" if err < 0.01 else "⚠️"
                st.caption(f"{cor} ref: {ref} · erro {err:.1e}")
            except (ValueError, ZeroDivisionError, TypeError):
                st.caption(f"ref: {ref} · {fonte}")


def mostrar_fonte(*objetos):
    """Aba Código: o trecho REAL de src/ que faz a conta."""
    for obj in objetos:
        st.caption(f"`{obj.__module__}.{obj.__name__}`")
        st.code(inspect.getsource(obj), language="python")


def carrega_csv(rel):
    import pandas as pd
    caminho = os.path.join(RAIZ, rel)
    if os.path.exists(caminho):
        return pd.read_csv(caminho)
    return None


def mostra_csv(rel, titulo=None):
    df = carrega_csv(rel)
    if df is None:
        st.info(f"{rel} ainda não foi gerado — rode `python main.py`.")
        return
    if titulo:
        st.caption(titulo)
    st.dataframe(df, width='stretch')

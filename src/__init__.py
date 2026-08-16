# -*- coding: utf-8 -*-
"""
Laboratório de espalhamento a baixa energia — IFSC/USP.

===============================================================================
 COMO OS ARQUIVOS ESTÃO ORGANIZADOS: POR NÚMERO DE CORPOS
===============================================================================

  src/comum/          o que serve para todo mundo
      constantes.py       constantes físicas e conversões de unidade
      solvers.py          integradores da equação radial (diferença central,
                          Numerov, grade alinhada à borda do poço)
      incerteza.py        extrapolação de Richardson e portões de validade

  src/um_corpo/       UMA partícula num potencial externo
      schrodinger.py      Schrödinger canônica: oscilador, hidrogênio, caixa

  src/dois_corpos/    DUAS partículas que interagem — o coração do laboratório
      potenciais.py       poço esférico, gaussiano, Pöschl-Teller, Lennard-Jones
      espalhamento.py     comprimento de espalhamento e alcance efetivo
      analitico.py        fórmulas fechadas para conferir o numérico
      ajuste.py           ajuste da expansão de alcance efetivo
      feshbach.py         ressonância de Feshbach: a(campo magnético)

  src/tres_corpos/    TRÊS partículas — física de Efimov
      trimero.py          equação hiper-radial, com e sem armadilha
      efimov.py           razões universais e o parâmetro de três corpos

  src/literatura/     comparação com valores publicados
      bench_referencias.py   potenciais e fórmulas de artigos específicos
      tabelas_artigo.py      as tabelas do artigo que reproduzimos

  src/relatorios/     geração de saída (figuras, relatórios, diário)
      gerar_figuras.py  gerar_relatorio.py  gerar_resultados.py
      gerar_inicio.py   diario.py

===============================================================================
 COMPATIBILIDADE
===============================================================================

Os dois jeitos de importar funcionam, e vão continuar funcionando:

    from src.tres_corpos import trimero     # explícito, diz onde mora
    from src import trimero                 # curto, o de sempre

O segundo funciona por causa dos apelidos registrados abaixo. Nenhum notebook,
teste ou script antigo precisa ser mudado.
"""

import importlib as _importlib
import sys as _sys

# nome curto -> subpacote onde o arquivo realmente está
ONDE_MORA = {
    "constantes": "comum",       "solvers": "comum",        "incerteza": "comum",
    "schrodinger": "um_corpo",
    "potenciais": "dois_corpos", "espalhamento": "dois_corpos",
    "analitico": "dois_corpos",  "ajuste": "dois_corpos",   "feshbach": "dois_corpos",
    "trimero": "tres_corpos",    "efimov": "tres_corpos",
    "bench_referencias": "literatura", "tabelas_artigo": "literatura",
    "gerar_figuras": "relatorios",  "gerar_inicio": "relatorios",
    "gerar_relatorio": "relatorios", "gerar_resultados": "relatorios",
    "diario": "relatorios",
}

# Os pesados (matplotlib, pandas) só carregam quando alguém pedir: __getattr__.
_LEVES = ("constantes", "solvers", "incerteza", "schrodinger", "potenciais",
          "espalhamento", "analitico", "ajuste", "feshbach", "trimero",
          "efimov", "tabelas_artigo", "bench_referencias")


def _apelidar(nome):
    """Carrega src.<subpacote>.<nome> e registra também como src.<nome>."""
    mod = _importlib.import_module(f".{ONDE_MORA[nome]}.{nome}", __name__)
    _sys.modules[f"{__name__}.{nome}"] = mod    # faz "from src.trimero import X" achar
    globals()[nome] = mod
    return mod


for _n in _LEVES:
    _apelidar(_n)


def __getattr__(nome):
    """PEP 562: 'from src import gerar_figuras' carrega só na hora."""
    if nome in ONDE_MORA:
        return _apelidar(nome)
    raise AttributeError(f"src não tem '{nome}'. Módulos: {sorted(ONDE_MORA)}")


__all__ = sorted(ONDE_MORA)

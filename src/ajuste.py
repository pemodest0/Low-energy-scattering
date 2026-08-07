# -*- coding: utf-8 -*-
"""
Ajuste dos parâmetros do potencial para reproduzir (a, r0) alvo
(Seção 4.5 do artigo): dois laços aninhados.

  laço interno : varia o parâmetro de INTENSIDADE (v, ou C6 no LJ)
                 até acertar o comprimento de espalhamento a;
  laço externo : varia o parâmetro de ALCANCE (mu, ou C12 no LJ)
                 até acertar o alcance efetivo r0.

Truque de robustez: perto da unitariedade a(v) diverge, mas 1/a(v) é
contínua e monótona ao cruzar o limiar de estado ligado.  Por isso o
laço interno resolve  1/a(v) = 1/a_alvo  (com 1/a_alvo = 0 na
unitariedade), usando scipy.optimize.brentq (v1; a v0 usava Illinois).

Observação (comentário de física): para os potenciais de um só fator de
escala (poço, mPT, gaussiano) vale a invariância de escala
a·mu = F(v), r0·mu = G(v), então a/r0 depende só de v e o laço externo
converge em poucas iterações.  Para o LJ não há essa invariância.

Ao final SEMPRE se verifica o número de nós (0 para a<0, 1 para o
dêuteron), como pede o artigo.
"""
import math
from dataclasses import dataclass

from . import espalhamento
from .potenciais import FABRICAS


# ------------------------------------------------------------ raiz 1D
# v1: usamos scipy.optimize.brentq (convergência superlinear garantida em
# intervalo com troca de sinal) no lugar do método de Illinois caseiro da
# v0.  O encapsulamento geométrico da raiz continua nosso.
from scipy.optimize import brentq


def _illinois(f, x0, x1, tol, max_iter=100):
    """Wrapper (nome herdado da v0): raiz de f em [x0, x1] via brentq.

    'tol' era tolerância em |f| na v0; brentq usa tolerância em x, então
    usamos xtol apertado e mantemos a interface.  Se f não troca de sinal
    no intervalo, brentq levanta ValueError (mesmo contrato da v0).
    """
    return brentq(f, x0, x1, xtol=1e-12, rtol=1e-14, maxiter=max_iter)


def _encapsula(f, x, fator=1.3, max_exp=40):
    """Expande geometricamente [x/fator^k, x*fator^m] até f trocar de sinal."""
    lo = hi = x
    flo = fhi = f(x)
    for _ in range(max_exp):
        if flo * fhi < 0:
            return (lo, hi) if lo < hi else (hi, lo)
        lo /= fator
        flo = f(lo)
        if flo * fhi < 0:
            return lo, hi
        hi *= fator
        fhi = f(hi)
    raise RuntimeError("não foi possível encapsular a raiz")


# ------------------------------------------------------------- ajuste
@dataclass
class ResultadoAjuste:
    p1: float           # intensidade (v ou C6)
    p2: float           # alcance (mu ou C12)
    a: float
    r0: float
    nos: int
    iteracoes: int
    convergiu: bool


def ajustar(nome_pot, a_alvo, r0_alvo, p1_ini, p2_ini,
            dr=2e-3, metodo="numerov",
            tol_a=1e-6, tol_r0=1e-4, max_iter_ext=30, nos_alvo=None):
    """Laços aninhados: interno resolve 1/a(p1) = 1/a_alvo; externo ajusta
    p2 até |r0 - r0_alvo| < tol_r0.  a_alvo = math.inf -> unitariedade."""
    fab = FABRICAS[nome_pot]
    inv_a_alvo = 0.0 if math.isinf(a_alvo) else 1.0 / a_alvo

    def resolve_p1(p2, p1_chute):
        """Laço interno: acha p1 tal que 1/a = 1/a_alvo (p2 fixo)."""
        f = lambda p1: 1.0 / espalhamento.calcular(fab(p1, p2), dr, metodo).a \
            - inv_a_alvo
        lo, hi = _encapsula(f, p1_chute)
        return _illinois(f, lo, hi, tol_a)

    p1, p2 = p1_ini, p2_ini
    it = 0
    for it in range(1, max_iter_ext + 1):
        p1 = resolve_p1(p2, p1)                       # acerta a
        res = espalhamento.calcular(fab(p1, p2), dr, metodo)
        if abs(res.r0 - r0_alvo) < tol_r0:
            break
        # laço externo: raiz de g(p2) = r0(p2; p1 re-ajustado) - r0_alvo
        def g(p2_teste):
            p1_loc = resolve_p1(p2_teste, p1)
            return espalhamento.calcular(fab(p1_loc, p2_teste),
                                         dr, metodo).r0 - r0_alvo
        lo, hi = _encapsula(g, p2)
        p2 = _illinois(g, lo, hi, tol_r0)
        p1 = resolve_p1(p2, p1)
        res = espalhamento.calcular(fab(p1, p2), dr, metodo)
        if abs(res.r0 - r0_alvo) < tol_r0:
            break

    convergiu = abs(res.r0 - r0_alvo) < tol_r0 and \
        abs(1.0 / res.a - inv_a_alvo) < 10 * tol_a
    # checagem final do número de nós (obrigatória — Seção 4.5)
    if nos_alvo is not None and res.nos != nos_alvo:
        convergiu = False
        print(f"  [aviso] nº de nós = {res.nos}, esperado {nos_alvo}!")
    return ResultadoAjuste(p1=p1, p2=p2, a=res.a, r0=res.r0,
                           nos=res.nos, iteracoes=it, convergiu=convergiu)

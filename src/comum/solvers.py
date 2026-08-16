# -*- coding: utf-8 -*-
"""
Integração da equação radial de energia zero.

FORMA UNIFICADA (v1.16) — uma equação para todas as dimensões e todos os
momentos angulares:

    u''(r) = [ 2 V(r) + (L^2 - 1/4)/r^2 ] u(r)      L = l + (d-2)/2

com a substituição u = r^{(d-1)/2} R.  Os quatro casos que interessam:

    d=3, onda parcial l  ->  L = l + 1/2  ->  L^2 - 1/4 = l(l+1)   <- o clássico
    d=2, momento m       ->  L = m        ->  L^2 - 1/4 = m^2 - 1/4
    d=1                  ->  L = l - 1/2  ->  0 (para l = 0)
    d=6, hiper-raio      ->  L = s        ->  L^2 - 1/4 = s^2 - 1/4  <- Efimov

A última linha é o ponto: o potencial efetivo -(s0^2 + 1/4)/R^2 da física
de Efimov NÃO é uma equação nova — é esta, com d = 6 e L = i s0.  O 1/4
vem do jacobiano em 6 dimensões (geometria), não da interação.

NOTA HONESTA: a *equação* é a mesma, mas o hiper-raio conserva solver
próprio (`src/efimov.py`), porque a torre de Efimov varre décadas de R e
exige a troca x = ln(R/R0) para caber numa grade uniforme.  A unificação
aqui cobre 1D/2D/3D e ondas parciais; o d=6 entra como verificação de
consistência (ver tests/test_solver_unificado.py).

Compatibilidade: o padrão é (l=0, d=3), onde L^2 - 1/4 = 0 exatamente e a
equação volta a ser u'' = 2 V u.  Nenhum resultado anterior muda.

Dois métodos (Seção 3.1 do artigo), escritos com o operador completo
W(r) = 2 V(r) + (L^2 - 1/4)/r^2:

  Método A - diferença central de 2a ordem (Eq. 99):
      u_{i+1} = 2 u_i - u_{i-1} + (dr)^2 W(r_i) u_i

  Método B - Numerov (Eq. 101), para y'' = -xi(x) y + s(x) com s = 0 e
  xi(r) = -W(r):
      u_{i+1} = [ 2 u_i (1 - 5 h2 xi_i) - u_{i-1} (1 + h2 xi_{i-1}) ]
                / (1 + h2 xi_{i+1}),    h2 = (dr)^2/12

Condições de contorno: u(r_min) = 0 e u(r_min + dr) = semente (a
normalização é livre, a equação é linear).  A semente segue a solução
regular r^{|L|+1/2}, o que reduz a contaminação pela solução irregular
r^{1/2-|L|} quando l > 0.  Integramos até R_match + dr, onde
R_match = r_min + N dr >= R é o primeiro ponto de grade além do alcance
(lá V ~ 0, então casar em R_match em vez de R não introduz erro).

Para o LJ o caroço repulsivo faz u crescer exponencialmente; usamos
re-escalonamento periódico (linearidade) para evitar overflow.
"""
import math
import numpy as np

ESCALA_MAX = 1e250   # limiar de re-escalonamento contra overflow


def grade(pot, dr):
    """Grade radial: r_i = r_min + i dr_eff, i = 0..N+1, com r_N = R exato.

    O passo pedido dr é levemente ajustado para dr_eff = (R - r_min)/N,
    de modo que a borda do potencial caia exatamente sobre um ponto da
    grade.  Isso elimina o erro de desalinhamento O(dr) na borda do poço
    esférico (sawtooth em função de mu) e não altera nada nos potenciais
    suaves.  O casamento com g0 = 1 - r/a é feito em r_N = R.
    """
    N = int(math.ceil((pot.R - pot.r_min) / dr - 1e-12))
    if N < 20:
        raise ValueError("Grade grossa demais: diminua dr.")
    dr_eff = (pot.R - pot.r_min) / N
    r = pot.r_min + dr_eff * np.arange(N + 2)
    return r, N, dr_eff


def L_efetivo(l=0, d=3):
    """L = l + (d-2)/2, para o momento angular USUAL em d dimensões.

    Vem de L^2 = lambda + (d-1)(d-3)/4 + 1/4 com o autovalor angular
    padrão lambda = l(l+d-2):

        L^2 = l(l+d-2) + [(d-1)(d-3)+1]/4 = (l + (d-2)/2)^2

    >>> L_efetivo(0, 3)**2 - 0.25        # onda-s em 3D
    0.0
    >>> round(L_efetivo(1, 3)**2 - 0.25, 12)   # onda-p: l(l+1) = 2
    2.0
    >>> L_efetivo(2, 2)                        # 2D, m = 2
    2.0

    ATENÇÃO — o caso hiperradial NÃO passa por aqui.  No problema de três
    corpos (d = 6) o autovalor angular NÃO é l(l+d-2): ele vem da equação
    de Faddeev/hiperangular e vale lambda = s^2 - 4.  Aí
    L^2 = s^2 - 4 + 15/4 + 1/4 = s^2, isto é, **L = s diretamente**.
    Use `integrar(..., L=s)` nesse caso, não `l=s, d=6`.
    """
    return l + (d - 2) / 2.0


def termo_centrifugo(r, L):
    """(L^2 - 1/4)/r^2, com r = 0 tratado como 0.

    Em r = 0 o termo diverge, mas ele sempre multiplica u(0) = 0; zerar
    ali evita 0*inf = nan sem alterar a solução.
    """
    c = L * L - 0.25
    if c == 0.0:
        return np.zeros_like(r)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(r > 0, c / np.where(r > 0, r, 1.0) ** 2, 0.0)
    return t


def integrar(pot, dr, metodo="numerov", l=0, d=3, L=None):
    """Integra u'' = [2V + (L^2-1/4)/r^2] u.

    L é obtido de (l, d) por `L_efetivo`, OU passado direto por `L=` — que
    é o caminho para o hiper-raio (L = s) e para qualquer caso em que o
    autovalor angular não seja o usual l(l+d-2).

    Padrão (l=0, d=3): o termo centrífugo é exatamente zero e a equação
    é a mesma de sempre, u'' = 2 V u.

    Devolve (r, u, N, V, dr_eff) com r[N] = R.  V é o potencial NU (sem o
    termo centrífugo), para que o cálculo de a e r0 não mude.
    """
    r, N, dr = grade(pot, dr)          # dr passa a ser o passo efetivo
    V = pot.V(r)
    L = L_efetivo(l, d) if L is None else float(L)
    W = 2.0 * V + termo_centrifugo(r, L)   # o operador completo

    u = np.empty(N + 2)
    u[0] = 0.0
    # semente: perto da origem as duas soluções vão como r^{1/2 +- L};
    # a regular (que anula em 0) é r^{|L|+1/2}  (= r^{l+1} em 3D).
    # Para potenciais com caroço (r_min > 0) qualquer semente serve.
    u[1] = (r[1] / r[-1]) ** (abs(L) + 0.5) if r[0] == 0.0 else 1.0
    if u[1] == 0.0:                        # underflow para L grande
        u[1] = 1e-300

    if metodo == "central":
        c = dr * dr
        for i in range(1, N + 1):
            u[i + 1] = 2.0 * u[i] - u[i - 1] + c * W[i] * u[i]
            if abs(u[i + 1]) > ESCALA_MAX:            # re-escalona (linearidade)
                u[: i + 2] /= ESCALA_MAX
    elif metodo == "numerov":
        xi = -W
        h2 = dr * dr / 12.0
        for i in range(1, N + 1):
            u[i + 1] = (2.0 * u[i] * (1.0 - 5.0 * h2 * xi[i])
                        - u[i - 1] * (1.0 + h2 * xi[i - 1])) \
                       / (1.0 + h2 * xi[i + 1])
            if abs(u[i + 1]) > ESCALA_MAX:
                u[: i + 2] /= ESCALA_MAX
    else:
        raise ValueError(f"método desconhecido: {metodo}")

    return r, u, N, V, dr

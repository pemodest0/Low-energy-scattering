# -*- coding: utf-8 -*-
"""
Os quatro potenciais de dois corpos do artigo (Seções 4.1-4.4), em unidades
adimensionais: hbar = m_r = 1, comprimentos em fm (l = 1 fm).

Nessas unidades o fator (hbar^2/m_r) = 1 e os potenciais ficam:

    Poço esférico : V(r) = -v * mu^2                 (r < R = 1/mu); 0 fora
    mPT           : V(r) = -v * mu^2 / cosh^2(mu r)
    Gaussiano     : V(r) = -v * mu^2 * exp(-r^2 mu^2)
    Lennard-Jones : ver nota de convenção na classe LennardJones

Alcance numérico R: menor r tal que |V(r)| <= EPS_CORTE = 1e-15
(para o poço, R é explícito).  Para o LJ também definimos r_min tal que
U(r_min) ~ 1e10 e impomos u(r) = 0 para r < r_min (caroço duro numérico).
"""
import math
import numpy as np

EPS_CORTE = 1e-15    # |V(R)| <= EPS_CORTE define o alcance numérico
V_CARROCO = 1e10     # U(r_min) ~ V_CARROCO define o corte do caroço do LJ


class Potencial:
    """Interface comum: V(r) vetorizado, alcance R, r_min e rótulos."""
    nome = "generico"
    r_min = 0.0        # início da integração (0, exceto para o LJ)

    def V(self, r):
        raise NotImplementedError

    def parametros(self):
        raise NotImplementedError


class PocoEsferico(Potencial):
    """V(r) = -v mu^2 para r < R = 1/mu; 0 para r >= R.  (Eq. 74 do artigo)"""
    nome = "poco"
    rotulo = "Poço esférico"

    def __init__(self, v, mu):
        self.v = float(v)
        self.mu = float(mu)
        self.R = 1.0 / self.mu           # alcance explícito

    def V(self, r):
        # O ponto de grade que cai exatamente na borda r = R é tratado como
        # DENTRO do poço (tolerância de ponto flutuante 1e-12 R): testes de
        # convergência mostram que essa convenção minimiza o erro O(dr)
        # associado à descontinuidade do potencial na borda.
        r = np.asarray(r, dtype=float)
        return np.where(r <= self.R * (1.0 + 1e-12), -self.v * self.mu**2, 0.0)

    def parametros(self):
        return {"v": self.v, "mu": self.mu}


class PoschlTeller(Potencial):
    """V(r) = -v mu^2 / cosh^2(mu r).  (Eq. 116 do artigo)"""
    nome = "mpt"
    rotulo = "Pöschl-Teller mod."

    def __init__(self, v, mu):
        self.v = float(v)
        self.mu = float(mu)
        # |V(R)| = eps  ->  cosh^2(mu R) = v mu^2/eps
        self.R = math.acosh(math.sqrt(self.v * self.mu**2 / EPS_CORTE)) / self.mu

    def V(self, r):
        r = np.asarray(r, dtype=float)
        return -self.v * self.mu**2 / np.cosh(self.mu * r)**2

    def parametros(self):
        return {"v": self.v, "mu": self.mu}


class Gaussiano(Potencial):
    """V(r) = -v mu^2 exp(-r^2 mu^2).  (Eq. 120 do artigo)"""
    nome = "gauss"
    rotulo = "Gaussiano"

    def __init__(self, v, mu):
        self.v = float(v)
        self.mu = float(mu)
        # |V(R)| = eps  ->  R = sqrt(ln(v mu^2/eps))/mu
        self.R = math.sqrt(math.log(self.v * self.mu**2 / EPS_CORTE)) / self.mu

    def V(self, r):
        r = np.asarray(r, dtype=float)
        return -self.v * self.mu**2 * np.exp(-(r * self.mu)**2)

    def parametros(self):
        return {"v": self.v, "mu": self.mu}


class LennardJones(Potencial):
    """Lennard-Jones (Eq. 121 do artigo), C12 em fm^10 e C6 em fm^4.

    NOTA IMPORTANTE DE CONVENÇÃO (verificada numericamente): a Eq. (121)
    do artigo escreve V = (hbar^2/m_r)[C12/r^12 - C6/r^6], mas as
    constantes tabeladas na Tabela 4 SÓ reproduzem os valores publicados
    de (a, r0) se a equação radial de energia zero for

        u''(r) = [C12/r^12 - C6/r^6] u(r),

    ou seja, com V = (hbar^2/(2 m_r))[C12/r^12 - C6/r^6].  Na nossa
    convenção u'' = 2 V u, isso significa V(r) = (1/2)(C12/r^12 - C6/r^6).
    (Provável deslize de convenção no artigo; ver relatório.)

    O potencial diverge em r -> 0: começamos a integrar em r_min, definido
    por U(r_min) = V_CARROCO, com u(r) = 0 para r < r_min (U = 2V).
    """
    nome = "lj"
    rotulo = "Lennard-Jones"

    def __init__(self, C12, C6):
        self.C12 = float(C12)
        self.C6 = float(C6)
        # cortes definidos sobre U(r) = C12/r^12 - C6/r^6 (adimensional):
        # r_min: U(r_min) = V_CARROCO (caroço repulsivo domina; Newton)
        self.r_min = self._resolve_U(V_CARROCO, (self.C12 / V_CARROCO)**(1.0 / 12))
        # R: |U(R)| = EPS_CORTE (cauda atrativa domina; Newton)
        self.R = self._resolve_U(-EPS_CORTE, (self.C6 / EPS_CORTE)**(1.0 / 6))

    def _resolve_U(self, alvo, chute):
        """Resolve U(r) = C12/r^12 - C6/r^6 = alvo por Newton."""
        r = chute
        for _ in range(100):
            f = self.C12 / r**12 - self.C6 / r**6 - alvo
            df = -12 * self.C12 / r**13 + 6 * self.C6 / r**7
            passo = f / df
            r -= passo
            if abs(passo) < 1e-14 * abs(r):
                break
        return r

    def V(self, r):
        r = np.asarray(r, dtype=float)
        rs = np.where(r > 0, r, np.inf)      # evita divisão por zero
        return 0.5 * (self.C12 / rs**12 - self.C6 / rs**6)

    def parametros(self):
        return {"C12": self.C12, "C6": self.C6}


# Fábricas usadas pelo módulo de ajuste (parâmetro1 = intensidade, 2 = alcance)
FABRICAS = {
    "poco": lambda p1, p2: PocoEsferico(p1, p2),
    "mpt": lambda p1, p2: PoschlTeller(p1, p2),
    "gauss": lambda p1, p2: Gaussiano(p1, p2),
    "lj": lambda p1, p2: LennardJones(p2, p1),   # (C6=intensidade, C12=alcance)
}

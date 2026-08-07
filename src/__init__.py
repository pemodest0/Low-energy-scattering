# -*- coding: utf-8 -*-
"""
Laboratório numérico de espalhamento de baixa energia (onda-s, energia zero).

Reprodução de: M. Macêdo-Lima e L. Madeira,
"Scattering length and effective range of microscopic two-body potentials",
Rev. Bras. Ens. Fís. 45, e20230079 (2023).
DOI: 10.1590/1806-9126-RBEF-2023-0079

Convenção de unidades (núcleo do programa):
    hbar = m_r = 1,  escala de comprimento l = 1 fm.
    epsilon = hbar^2 / (m_r l^2)  -> unidade de energia.
    Equação resolvida:  u''(rbar) = 2 Vbar(rbar) u(rbar)   (energia zero, l=0)
"""

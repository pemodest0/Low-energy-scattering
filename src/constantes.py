# -*- coding: utf-8 -*-
"""
Constantes físicas usadas na validação da Tabela 1 do artigo
(energias de estado ligado do dímero de 4He e do dêuteron).

O núcleo numérico do laboratório é adimensional (hbar = m_r = 1, l = 1 fm);
estas constantes só entram na hora de "recuperar as unidades".
"""
# ---------------------------------------------------------------- CODATA
HBARC_MEV_FM = 197.3269804          # hbar*c em MeV.fm
MP_MEV = 938.27208816               # massa do próton, MeV/c^2
MN_MEV = 939.56542052               # massa do nêutron, MeV/c^2

HBAR_SI = 1.054571817e-34           # J.s
U_KG = 1.66053906660e-27            # unidade de massa atômica, kg
KB_SI = 1.380649e-23                # J/K
M_HE4_U = 4.002602                  # massa do 4He em u

# Massas reduzidas -------------------------------------------------------
M_R_NP_MEV = MP_MEV * MN_MEV / (MP_MEV + MN_MEV)   # dêuteron (n-p)
M_R_NN_MEV = MN_MEV / 2.0                          # par nêutron-nêutron
M_R_HE_U = M_HE4_U / 2.0                           # dímero de 4He (em u)

# hbar^2/(2 m_r) em unidades convenientes --------------------------------
def hbar2_over_2mr_MeV_fm2(m_r_mev):
    """hbar^2/(2 m_r) em MeV.fm^2, com m_r em MeV/c^2."""
    return HBARC_MEV_FM**2 / (2.0 * m_r_mev)

def hbar2_over_2mr_K_A2(m_r_u):
    """hbar^2/(2 m_r) em K.Angstrom^2 (dividido por k_B), com m_r em u."""
    joule_A2 = HBAR_SI**2 / (2.0 * m_r_u * U_KG * (1e-10)**2)
    return joule_A2 / KB_SI

EULER_GAMMA = 0.5772156649015328606  # constante de Euler-Mascheroni

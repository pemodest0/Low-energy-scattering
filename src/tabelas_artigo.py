# -*- coding: utf-8 -*-
"""
Valores-alvo e de referência das Tabelas 1-4 do artigo (unidades: fm;
C12 em fm^10, C6 em fm^4).  Usados para validação e como chutes iniciais
do ajuste (Seção 4.5).
"""
import math

# ------------------------------------------------ Tabela 2: alvos físicos
TABELA2 = {
    "nn":         {"a": -18.5,     "r0": 2.7, "nos": 0},
    "unitario":   {"a": math.inf,  "r0": 1.0, "nos": 0},   # |a| -> infinito
    "deuteron":   {"a": 5.4,       "r0": 1.7, "nos": 1},
}
ROTULOS_CASO = {"nn": "Nêutron-nêutron", "unitario": "Unitariedade",
                "deuteron": "Dêuteron"}

# --------------------------- Tabela 3: parâmetros publicados (poço/mPT/gauss)
# (p1, p2) = (v, mu); 'a_ref' e 'r0_ref' são os valores reportados no artigo.
TABELA3 = {
    ("nn", "poco"):        {"p1": 1.1096, "p2": 0.3918, "a_ref": -18.52, "r0_ref": 2.70},
    ("nn", "mpt"):         {"p1": 0.9071, "p2": 0.7991, "a_ref": -18.51, "r0_ref": 2.70},
    ("nn", "gauss"):       {"p1": 1.2121, "p2": 0.5672, "a_ref": -18.55, "r0_ref": 2.70},
    ("unitario", "poco"):  {"p1": 1.2337, "p2": 1.0000, "a_ref": -1e5,   "r0_ref": 1.00},
    ("unitario", "mpt"):   {"p1": 1.0000, "p2": 2.0000, "a_ref": 1e9,    "r0_ref": 1.00},
    ("unitario", "gauss"): {"p1": 1.3420, "p2": 1.4349, "a_ref": -1e5,   "r0_ref": 1.00},
    ("deuteron", "poco"):  {"p1": 1.7575, "p2": 0.5000, "a_ref": 5.40,   "r0_ref": 1.70},
    ("deuteron", "mpt"):   {"p1": 1.4388, "p2": 0.8631, "a_ref": 5.40,   "r0_ref": 1.73},
    ("deuteron", "gauss"): {"p1": 1.9102, "p2": 0.6754, "a_ref": 5.40,   "r0_ref": 1.70},
}

# ----------------------------------- Tabela 4: Lennard-Jones (C6, C12)
# (p1, p2) = (C6, C12) na convenção do módulo de ajuste.
TABELA4 = {
    ("nn", "lj"):       {"p1": 9.86668911, "p2": 3.08836698, "a_ref": -18.5, "r0_ref": 2.71},
    ("unitario", "lj"): {"p1": 0.26462461, "p2": 0.00034068, "a_ref": -1e5,  "r0_ref": 1.00},
    ("deuteron", "lj"): {"p1": 6.81472000, "p2": 0.90485319, "a_ref": 5.4,   "r0_ref": 1.70},
}

# --------------------------- Tabela 1: energias de estado ligado (validação)
TABELA1 = {
    "he4_dimer": {"a": 90.4, "r0": 8.0,           # Angstrom
                  "E_ref": -1.62e-3,              # K
                  "E_zr_ref": -1.48e-3, "E_fr_ref": -1.63e-3},
    "deuteron":  {"a": 5.4112, "r0": 1.7436,      # fm
                  "E_ref": -2.224,                # MeV
                  "E_zr_ref": -1.416, "E_fr_ref": -2.223},
}

NOMES_PARAM = {"poco": ("v", "mu"), "mpt": ("v", "mu"),
               "gauss": ("v", "mu"), "lj": ("C6", "C12")}

# -*- coding: utf-8 -*-
"""
Gera os resultados numéricos do laboratório e escreve CSVs, fragmentos
LaTeX e um resumo em resultados/.

Etapas:
  1. validação analítica (poço e mPT) - dois métodos;
  2. reprodução das Tabelas 3 e 4 a partir dos parâmetros publicados;
  3. energias de estado ligado da Tabela 1;
  4. ajuste dos parâmetros (Seção 4.5) para todos os casos;
  5. estudo de convergência em dr (para a figura de convergência).
"""
import csv
import math
import os

from . import espalhamento, analitico, ajuste
from .potenciais import FABRICAS
from .tabelas_artigo import (TABELA1, TABELA2, TABELA3, TABELA4,
                             NOMES_PARAM)
from . import constantes as cte

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_RES = os.path.join(AQUI, "resultados")


def _erro_rel(x, ref):
    return abs(x - ref) / abs(ref) if ref not in (0, None) else float("nan")


# ------------------------------------------------- 1. validação analítica
def validacao_analitica(dr=1e-3):
    """Poço e mPT: numérico (2 métodos) vs fórmulas fechadas."""
    linhas = []
    casos = [
        ("poco", "nn"), ("poco", "unitario"), ("poco", "deuteron"),
        ("mpt", "nn"), ("mpt", "unitario"), ("mpt", "deuteron"),
    ]
    for pot, caso in casos:
        p = TABELA3[(caso, pot)]
        v, mu = p["p1"], p["p2"]
        if pot == "poco":
            R = 1.0 / mu
            a_ana = analitico.a_poco(v, R)
            r0_ana = analitico.r0_poco(v, R)
        else:
            a_ana = analitico.a_mpt(v, mu)
            r0_ana = analitico.r0_mpt_unitario(mu) if caso == "unitario" else None
        for met in ("central", "numerov"):
            res = espalhamento.calcular(FABRICAS[pot](v, mu), dr=dr, metodo=met)
            linhas.append({
                "potencial": pot, "caso": caso, "metodo": met, "dr": dr,
                "a_num": res.a, "a_analitico": a_ana,
                "erro_rel_a": _erro_rel(res.a, a_ana) if abs(a_ana) < 1e4 else abs(1/res.a - 1/a_ana),
                "r0_trap": res.r0_trap, "r0_simpson": res.r0_simp,
                "r0_analitico": r0_ana,
                "erro_rel_r0": _erro_rel(res.r0_simp, r0_ana) if r0_ana else None,
                "nos": res.nos,
            })
    return linhas


# ------------------------------------- 2. reprodução das Tabelas 3 e 4
def reproducao_tabelas(dr=1e-3, metodo="numerov"):
    """Roda os parâmetros PUBLICADOS e compara (a, r0) com o artigo."""
    linhas = []
    for tabela, nome_tab in ((TABELA3, "Tabela 3"), (TABELA4, "Tabela 4")):
        for (caso, pot), p in tabela.items():
            fab = FABRICAS[pot]
            # convenção da fábrica: p1 = intensidade (v ou C6), p2 = alcance
            if pot == "lj":
                res = espalhamento.calcular(fab(p["p1"], p["p2"]), dr, metodo)
            else:
                res = espalhamento.calcular(fab(p["p1"], p["p2"]), dr, metodo)
            alvo = TABELA2[caso]
            unit = math.isinf(alvo["a"])
            linhas.append({
                "tabela": nome_tab, "caso": caso, "potencial": pot,
                NOMES_PARAM[pot][0]: p["p1"], NOMES_PARAM[pot][1]: p["p2"],
                "a_num": res.a, "a_artigo": p["a_ref"],
                "erro_rel_a": abs(1.0 / res.a) if unit else _erro_rel(res.a, alvo["a"]),
                "r0_num": res.r0_simp, "r0_artigo": p["r0_ref"],
                "erro_rel_r0": _erro_rel(res.r0_simp, alvo["r0"]),
                "nos": res.nos, "nos_esperado": alvo["nos"],
            })
    return linhas


# --------------------------------------------- 3. energias da Tabela 1
def energias_tabela1():
    linhas = []
    # dímero de 4He (Angstrom, Kelvin)
    d = TABELA1["he4_dimer"]
    h22m = cte.hbar2_over_2mr_K_A2(cte.M_R_HE_U)
    linhas.append({
        "sistema": "dimero_4He", "a": d["a"], "r0": d["r0"],
        "unid": "Angstrom/K", "E_ref": d["E_ref"],
        "E_zr": analitico.energia_zr(d["a"], h22m),
        "E_zr_artigo": d["E_zr_ref"],
        "E_fr": analitico.energia_fr(d["a"], d["r0"], h22m),
        "E_fr_artigo": d["E_fr_ref"],
    })
    # dêuteron (fm, MeV)
    d = TABELA1["deuteron"]
    h22m = cte.hbar2_over_2mr_MeV_fm2(cte.M_R_NP_MEV)
    linhas.append({
        "sistema": "deuteron", "a": d["a"], "r0": d["r0"],
        "unid": "fm/MeV", "E_ref": d["E_ref"],
        "E_zr": analitico.energia_zr(d["a"], h22m),
        "E_zr_artigo": d["E_zr_ref"],
        "E_fr": analitico.energia_fr(d["a"], d["r0"], h22m),
        "E_fr_artigo": d["E_fr_ref"],
    })
    return linhas


# ------------------------------------------------- 4. ajuste (Seção 4.5)
# chutes iniciais genéricos (deslocados dos valores publicados de propósito,
# para demonstrar a convergência do algoritmo de laços aninhados)
def rodar_ajustes(dr=2e-3, verbose=True):
    linhas = []
    for tabela in (TABELA3, TABELA4):
        for (caso, pot), p in tabela.items():
            alvo = TABELA2[caso]
            p1_ini, p2_ini = 1.15 * p["p1"], 0.85 * p["p2"]
            if verbose:
                print(f"ajustando {pot}/{caso} ...", flush=True)
            try:
                aj = ajuste.ajustar(pot, alvo["a"], alvo["r0"],
                                    p1_ini, p2_ini, dr=dr,
                                    nos_alvo=alvo["nos"])
                n1, n2 = NOMES_PARAM[pot]
                linhas.append({
                    "caso": caso, "potencial": pot,
                    f"{n1}_ajustado": aj.p1, f"{n1}_artigo": p["p1"],
                    f"erro_rel_{n1}": _erro_rel(aj.p1, p["p1"]),
                    f"{n2}_ajustado": aj.p2, f"{n2}_artigo": p["p2"],
                    f"erro_rel_{n2}": _erro_rel(aj.p2, p["p2"]),
                    "a_obtido": aj.a, "r0_obtido": aj.r0,
                    "nos": aj.nos, "convergiu": aj.convergiu,
                })
            except Exception as exc:       # v0: não deixar um caso derrubar tudo
                linhas.append({"caso": caso, "potencial": pot,
                               "convergiu": False, "erro": str(exc)})
    return linhas


# ------------------------------------------------- 5. convergência em dr
def convergencia():
    """Erro relativo de a e r0 vs dr, comparado com resultados analíticos.
    a: poço (nn) e mPT (nn);  r0: poço (nn) e mPT unitário (r0 = 2/mu)."""
    drs = [2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4]
    linhas = []
    cfgs = [
        ("poco", 1.1096, 0.3918, "nn"),
        ("mpt", 0.9071, 0.7991, "nn"),
        ("mpt", 1.0, 2.0, "unitario"),
    ]
    for pot, v, mu, caso in cfgs:
        if pot == "poco":
            a_ana = analitico.a_poco(v, 1.0 / mu)
            r0_ana = analitico.r0_poco(v, 1.0 / mu)
        else:
            a_ana = analitico.a_mpt(v, mu)
            r0_ana = analitico.r0_mpt_unitario(mu) if caso == "unitario" else None
        for met in ("central", "numerov"):
            for dr in drs:
                res = espalhamento.calcular(FABRICAS[pot](v, mu), dr, met)
                linhas.append({
                    "potencial": pot, "caso": caso, "metodo": met, "dr": dr,
                    "N": len(res.r) - 2,
                    "a": res.a,
                    "erro_a": (_erro_rel(res.a, a_ana)
                               if abs(a_ana) < 1e6 else abs(1/res.a - 1/a_ana)),
                    "r0_trap": res.r0_trap, "r0_simpson": res.r0_simp,
                    "erro_r0_trap": _erro_rel(res.r0_trap, r0_ana) if r0_ana else None,
                    "erro_r0_simpson": _erro_rel(res.r0_simp, r0_ana) if r0_ana else None,
                })
    return linhas


# ----------------------------------------------------------- escrita
def _escreve_csv(nome, linhas):
    os.makedirs(DIR_RES, exist_ok=True)
    caminho = os.path.join(DIR_RES, nome)
    chaves = []
    for ln in linhas:                       # união das chaves, ordem estável
        for k in ln:
            if k not in chaves:
                chaves.append(k)
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=chaves)
        w.writeheader()
        w.writerows(linhas)
    return caminho


def gerar_tudo(dr=1e-3, dr_ajuste=2e-3):
    print("[1/5] validação analítica ...", flush=True)
    va = validacao_analitica(dr)
    _escreve_csv("validacao_analitica.csv", va)

    print("[2/5] reprodução das Tabelas 3 e 4 ...", flush=True)
    rt = reproducao_tabelas(dr)
    _escreve_csv("reproducao_tabelas_3_4.csv", rt)

    print("[3/5] energias da Tabela 1 ...", flush=True)
    e1 = energias_tabela1()
    _escreve_csv("energias_tabela1.csv", e1)

    print("[4/5] ajustes (Seção 4.5) ...", flush=True)
    aj = rodar_ajustes(dr_ajuste)
    _escreve_csv("ajustes.csv", aj)

    print("[5/5] convergência ...", flush=True)
    cv = convergencia()
    _escreve_csv("convergencia.csv", cv)

    _escreve_resumo(va, rt, e1, aj)
    return va, rt, e1, aj, cv


def _escreve_resumo(va, rt, e1, aj):
    """resumo.md legível, com destaque para os erros relativos."""
    L = []
    L.append("# Resumo dos resultados — laboratório de espalhamento\n")
    L.append("Reprodução de Macêdo-Lima & Madeira, RBEF 45, e20230079 (2023).\n")

    L.append("\n## Validação analítica (poço esférico e mPT)\n")
    L.append("| potencial | caso | método | a_num (fm) | a_analítico | erro rel a | r0_Simpson | r0_analítico | erro rel r0 | nós |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for x in va:
        r0a = f"{x['r0_analitico']:.6f}" if x["r0_analitico"] else "—"
        er0 = f"{x['erro_rel_r0']:.2e}" if x["erro_rel_r0"] is not None else "—"
        L.append(f"| {x['potencial']} | {x['caso']} | {x['metodo']} | "
                 f"{x['a_num']:.5g} | {x['a_analitico']:.5g} | {x['erro_rel_a']:.2e} | "
                 f"{x['r0_simpson']:.6f} | {r0a} | {er0} | {x['nos']} |")

    L.append("\n## Reprodução das Tabelas 3 e 4 (parâmetros publicados)\n")
    L.append("| tabela | caso | potencial | a_num | a_alvo | erro rel a | r0_num | r0_alvo | erro rel r0 | nós (esperado) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for x in rt:
        alvo = TABELA2[x["caso"]]
        a_alvo = "±inf" if math.isinf(alvo["a"]) else f"{alvo['a']}"
        L.append(f"| {x['tabela']} | {x['caso']} | {x['potencial']} | "
                 f"{x['a_num']:.4g} | {a_alvo} | {x['erro_rel_a']:.2e} | "
                 f"{x['r0_num']:.4f} | {alvo['r0']} | {x['erro_rel_r0']:.2e} | "
                 f"{x['nos']} ({x['nos_esperado']}) |")

    L.append("\n## Energias de estado ligado (Tabela 1)\n")
    L.append("| sistema | E_zr | E_zr artigo | E_fr | E_fr artigo | E exp. |")
    L.append("|---|---|---|---|---|---|")
    for x in e1:
        L.append(f"| {x['sistema']} | {x['E_zr']:.4g} | {x['E_zr_artigo']} | "
                 f"{x['E_fr']:.4g} | {x['E_fr_artigo']} | {x['E_ref']} |")

    L.append("\n## Ajuste de parâmetros (Seção 4.5) vs artigo\n")
    L.append("| caso | potencial | p1 ajustado | p1 artigo | erro rel | p2 ajustado | p2 artigo | erro rel | a obtido | r0 obtido | nós | convergiu |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for x in aj:
        if not x.get("convergiu", False) and "erro" in x:
            L.append(f"| {x['caso']} | {x['potencial']} | FALHOU: {x['erro']} |")
            continue
        n1, n2 = NOMES_PARAM[x["potencial"]]
        L.append(f"| {x['caso']} | {x['potencial']} | "
                 f"{x[f'{n1}_ajustado']:.5f} | {x[f'{n1}_artigo']} | {x[f'erro_rel_{n1}']:.2e} | "
                 f"{x[f'{n2}_ajustado']:.5f} | {x[f'{n2}_artigo']} | {x[f'erro_rel_{n2}']:.2e} | "
                 f"{x['a_obtido']:.4g} | {x['r0_obtido']:.4f} | {x['nos']} | {x['convergiu']} |")

    L.append("\n> Nota de convenção: os C12/C6 da Tabela 4 do artigo reproduzem os alvos "
             "somente com u'' = (C12/r^12 - C6/r^6) u, isto é, V_LJ = (hbar^2/2m_r)[...] "
             "e não (hbar^2/m_r)[...] como impresso na Eq. (121). Ver src/potenciais.py.\n")

    os.makedirs(DIR_RES, exist_ok=True)
    with open(os.path.join(DIR_RES, "resumo.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    gerar_tudo()

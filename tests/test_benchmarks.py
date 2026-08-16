# -*- coding: utf-8 -*-
"""
Testes dirigidos por dados: cada entrada de `referencias/BENCHMARKS.yaml`
que tenha um bloco `teste:` vira um caso de teste automático.

Filosofia: um valor de literatura sem teste é folclore com DOI. Se o
código deixar de reproduzir um número publicado, o pytest grita.

Rode:   python -m pytest tests/test_benchmarks.py -v
Só um:  python -m pytest tests/test_benchmarks.py -k efimov -v
"""
import math
import pathlib

import pytest
import yaml

from src import espalhamento
from src import constantes as cte
from src.potenciais import PocoEsferico, PoschlTeller, Gaussiano, LennardJones

RAIZ = pathlib.Path(__file__).resolve().parents[1]
ARQ = RAIZ / "referencias" / "BENCHMARKS.yaml"

FABRICAS = {
    "poco": PocoEsferico,
    "mpt": PoschlTeller,
    "gauss": Gaussiano,
    "lj": LennardJones,
}


def _carregar():
    with open(ARQ, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


DADOS = _carregar()
COM_TESTE = [e for e in DADOS["entradas"] if e.get("teste")]
SEM_TESTE = [e for e in DADOS["entradas"] if not e.get("teste")]


# ------------------------------------------------------ integridade do YAML
def test_todo_valor_tem_convencao():
    """Regra do repositório: nenhum número sem convenção declarada."""
    faltando = [e["id"] for e in DADOS["entradas"] if not e.get("convencao")]
    assert not faltando, f"entradas sem convenção: {faltando}"


def test_toda_fonte_existe():
    conhecidas = set(DADOS["fontes"]) | {"PENDENTE"}
    orfas = [e["id"] for e in DADOS["entradas"] if e.get("fonte") not in conhecidas]
    assert not orfas, f"entradas com fonte desconhecida: {orfas}"


def test_ids_unicos():
    ids = [e["id"] for e in DADOS["entradas"]]
    assert len(ids) == len(set(ids))


# ------------------------------------------------------ auxiliares de física
def _m_r(rotulo):
    return {"np": cte.M_R_NP_MEV, "nn": cte.M_R_NN_MEV, "he": cte.M_R_HE_U}[rotulo]


def _pref(rotulo):
    """hbar^2/(2 m_r) na unidade certa: MeV.fm^2 (nuclear) ou K.A^2 (He)."""
    if rotulo == "he":
        return cte.hbar2_over_2mr_K_A2(_m_r(rotulo))
    return cte.hbar2_over_2mr_MeV_fm2(_m_r(rotulo))


def _escala(rotulo, x):
    """He: converte K -> mK. Nuclear: já está em MeV."""
    return x * 1e3 if rotulo == "he" else x


def energia_zr(a, rot):
    """Eq. (93): E_zr = -hbar^2/(2 m_r a^2)."""
    return _escala(rot, -_pref(rot) / a**2)


def energia_fr(a, r0, rot):
    """Eq. (95): raiz de 1/a = kappa - r0 kappa^2/2."""
    kappa = (1.0 - math.sqrt(1.0 - 2.0 * r0 / a)) / r0
    return _escala(rot, -_pref(rot) * kappa**2)


# ------------------------------------------------------ o teste parametrizado
def _param(e):
    """Entrada marcada com `esperado: diverge` vira xfail — divergência é
    resultado registrado, não teste quebrado."""
    if e.get("esperado") == "diverge":
        return pytest.param(e, marks=pytest.mark.xfail(
            reason=f"divergência registrada: {e['id']} (ver CONVENCOES.md)",
            strict=False))
    return e


@pytest.mark.parametrize("entrada", [_param(e) for e in COM_TESTE],
                         ids=[e["id"] for e in COM_TESTE])
def test_benchmark(entrada):
    t = entrada["teste"]
    tipo = t["tipo"]
    alvo = entrada["valor"]
    rtol = float(t.get("rtol", 1e-3))

    if tipo == "espalhamento":
        pot = FABRICAS[t["potencial"]](*t["params"])
        res = espalhamento.calcular(pot, dr=1e-3)
        obtido = {"a": res.a, "r0": res.r0}[t["obs"]]

    elif tipo == "energia_zr":
        obtido = energia_zr(t["a"], t["m_r"])

    elif tipo == "energia_fr":
        obtido = energia_fr(t["a"], t["r0"], t["m_r"])

    elif tipo == "efimov_s0":
        from src import efimov
        obtido = efimov.s0_universal()

    elif tipo == "efimov_razao_E":
        from src import efimov
        obtido = math.exp(2.0 * math.pi / efimov.s0_universal())

    elif tipo == "feshbach_B0":
        from src import feshbach
        obtido = feshbach.RESSONANCIAS_39K[t["estado"]][t["indice"]]["B0"]

    elif tipo == "feshbach_lvdw":
        from src import feshbach
        obtido = feshbach.comprimento_vdw_a0()

    else:
        pytest.fail(f"tipo de teste desconhecido: {tipo}")

    assert obtido == pytest.approx(alvo, rel=rtol), (
        f"\n  benchmark : {entrada['id']}"
        f"\n  publicado : {alvo} {entrada.get('unidade','')}"
        f"\n  obtido    : {obtido}"
        f"\n  convenção : {entrada['convencao']}  (ver referencias/CONVENCOES.md)"
    )


# ------------------------------------------------------ inventário (não falha)
def test_inventario(capsys):
    """Não falha nunca: só imprime o que ainda falta. Use -s para ver."""
    pend = [e for e in DADOS["entradas"] if e.get("acao_humana")]
    with capsys.disabled():
        print(f"\n  BENCHMARKS: {len(DADOS['entradas'])} entradas | "
              f"{len(COM_TESTE)} testadas | {len(SEM_TESTE)} sem teste | "
              f"{len(pend)} aguardando ação humana")
        for e in pend:
            print(f"    [ ] {e['id']:<28} {e['acao_humana']}")

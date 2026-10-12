# -*- coding: utf-8 -*-
"""Fase 0 da segunda revisao: medir tudo o que esta em disputa. Nao altera
nada, so imprime.

Rodar de dentro de notebooks/2_dois_corpos:

    python fase0_medir.py

Cada bloco imprime o valor do PDF, o valor medido, e o veredito. Onde o Lucas
deu um numero, ele aparece na coluna de comparacao.
"""
import math
import sys
import numpy as np

import lab


LARGURA = 78


def titulo(n, texto):
    print("\n" + "=" * LARGURA)
    print(f"  {n}. {texto}")
    print("=" * LARGURA)


def veredito(nosso, deles, tol=0.02):
    if deles is None:
        return ""
    return "  <- LUCAS TEM RAZAO" if abs(nosso / deles - 1.0) > tol else "  (bate)"


# ===========================================================================
titulo(1, "A amplificacao: PDF diz 48.7, Lucas mede 32.6")
# ===========================================================================
# d ln a / d ln(hbar^2/2mu) para o HFD-B, por diferenca central em varias
# larguras de passo. Se o valor nao estabilizar, o problema e a derivada.
h2 = lab.H2_2MU_HE4
a_de = lambda h: lab.scattering(lab.AzizHFDB(h))[0]
a0 = a_de(h2)
print(f"  h2_2mu = {h2}   ->   a = {a0:.6f}")
print("\n    passo relativo    derivada logaritmica")
derivadas = []
for d in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
    hp, hm = h2 * (1.0 + d), h2 * (1.0 - d)
    der = (math.log(a_de(hp)) - math.log(a_de(hm))) / (math.log(hp) - math.log(hm))
    derivadas.append(der)
    print(f"      {d:.0e}            {der:8.3f}")
melhor = derivadas[-1]
print(f"\n    PDF diz            48.7")
print(f"    Lucas mede         32.6  (e 31.5 pelos dois valores finitos do texto)")
print(f"    medido agora     {melhor:6.1f}{veredito(melhor, 48.7, 0.05)}")

# a "razao efetiva pelos dois valores finitos" que ele cita: reproduzir tambem
h_cod = 12.11932
h_pub = 12.12
if abs(h_cod - h_pub) > 0:
    r = (math.log(a_de(h_pub)) - math.log(a_de(h_cod))) / (math.log(h_pub) - math.log(h_cod))
    print(f"    razao pelos dois valores do texto: {r:.1f}")


# ===========================================================================
titulo(2, "Apendice B: a tabela de tuning bate com o notebook?")
# ===========================================================================
# Lucas: nn/lj deve dar C6=9.768407, C12=3.011567 (PDF diz 9.760657, 3.005560)
#        deuteron/lj deve dar C6=6.883344, C12=0.926958 (PDF: 6.846862, 0.915045)
print("  sistema    potencial     p1 (C6 no lj)     p2 (C12 no lj)    ok?")
saida = {}
for sistema, alvo in lab.TARGETS.items():
    for nome in ("well", "mpt", "gauss", "lj"):
        pub = lab.PUBLISHED.get((sistema, nome))
        if pub is None:
            continue
        try:
            p1, p2, a, r0, nodes, ok = lab.tune(
                nome, alvo["a"], alvo["r0"], pub["p1"], pub["p2"],
                nodes_target=alvo["nodes"])
            saida[(sistema, nome)] = (p1, p2)
            print(f"  {sistema:10s} {nome:10s} {p1:14.6f} {p2:16.6f}    "
                  f"{'sim' if ok else 'NAO'}")
        except Exception as e:
            print(f"  {sistema:10s} {nome:10s}  FALHOU: {type(e).__name__}: {e}")

print("\n  confronto com os numeros do Lucas:")
for chave, (c6L, c12L) in ((("nn", "lj"), (9.768407, 3.011567)),
                           (("deuteron", "lj"), (6.883344, 0.926958))):
    if chave in saida:
        c6, c12 = saida[chave]
        print(f"    {chave[0]:10s} C6  medido {c6:10.6f}  Lucas {c6L:10.6f}"
              f"{veredito(c6, c6L, 1e-5)}")
        print(f"    {chave[0]:10s} C12 medido {c12:10.6f}  Lucas {c12L:10.6f}"
              f"{veredito(c12, c12L, 1e-5)}")


# ===========================================================================
titulo(3, "Lennard-Jones: r_min usa so o termo repulsivo. Quanto isso custa?")
# ===========================================================================
# O codigo resolve 0.5 C12 / r^10 = EPS_CORE, ignorando -C6/r^6.
# Medir |V(r_min)| r_min^2 de verdade, com o potencial completo.
from scipy.optimize import brentq
print("  sistema     |V(r_min)| r_min^2    EPS_CORE     desvio")
for sistema in ("nn", "unitarity", "deuteron"):
    pub = lab.PUBLISHED.get((sistema, "lj"))
    if pub is None:
        continue
    p = lab.LennardJones(pub["p1"], pub["p2"])
    real = abs(p.V(p.r_min)) * p.r_min ** 2
    print(f"  {sistema:10s} {real:18.4f} {lab.EPS_CORE:12.0f} "
          f"{(real/lab.EPS_CORE - 1.0)*100:+8.2f}%")

    # e onde estaria r_min se resolvido com o potencial completo
    f = lambda r: abs(p.V(r)) * r * r - lab.EPS_CORE
    exato = brentq(f, p.r_min * 0.5, p.r_min * 1.5, xtol=1e-15)
    print(f"             r_min fechado {p.r_min:.8f}   completo {exato:.8f}   "
          f"desloca {(exato/p.r_min - 1.0)*100:+.3f}%")


# ===========================================================================
titulo(4, "Lennard-Jones para helio-4: qual a energia? (PDF deixa em branco)")
# ===========================================================================
# Lucas reporta -1.627 mK. Reproduzir pelo mesmo caminho do notebook.
try:
    escala = None
    hfdb = lab.AzizHFDB(lab.H2_2MU_HE4)
    a_he, r0_he, _ = lab.scattering(hfdb)
    print(f"  alvo do helio pelo HFD-B: a = {a_he:.4f}, r0 = {r0_he:.4f}")
    p1, p2 = lab.GUESS["lj"]
    s, sc, a, r0, nodes, ok = lab.tune("lj", a_he, r0_he, p1, p2, nodes_target=1)
    pot = lab.LennardJones(s, sc)
    guess = lab.E_finite_range(a, r0, lab.H2_2MU_HE4) / (2 * lab.H2_2MU_HE4)
    E_mK = lab.bound_energy(pot, guess) * 2 * lab.H2_2MU_HE4 * 1e3
    print(f"  sintonizado: C6 = {s:.6f}, C12 = {sc:.6f}, nos = {nodes}, ok = {ok}")
    print(f"  energia medida  {E_mK:.4f} mK")
    print(f"  Lucas reporta   -1.627 mK{veredito(E_mK, -1.627, 0.01)}")
except Exception as e:
    print(f"  FALHOU: {type(e).__name__}: {e}")


# ===========================================================================
titulo(5, "Secao 15 contra Secao 19: EPS_TAIL = 1e-12 aproxima ou desestabiliza?")
# ===========================================================================
# A tabela da Sec 15 diz que 1e-12 aproxima o resultado a 0.004%.
# A Sec 19 diz que ali ja esta instavel. So uma pode estar certa.
import importlib
base = {}
print("  EPS_TAIL      poco: a, r0            Lennard-Jones: a, r0")
for eps in (1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14):
    lab.EPS_TAIL = eps
    linha = f"  {eps:.0e}   "
    for sistema, nome in (("deuteron", "well"), ("deuteron", "lj")):
        pub = lab.PUBLISHED[(sistema, nome)]
        p = lab.POTENTIALS[nome](pub["p1"], pub["p2"])
        a, r0, _ = lab.scattering(p)
        linha += f"  {a:10.6f} {r0:10.6f} "
        if eps == 1e-6:
            base[(sistema, nome)] = (a, r0)
    print(linha)

lab.EPS_TAIL = 1e-6
print("\n  desvio relativo em r0 contra EPS_TAIL = 1e-6 (o valor atual):")
for eps in (1e-10, 1e-12, 1e-14):
    lab.EPS_TAIL = eps
    for sistema, nome in (("deuteron", "well"), ("deuteron", "lj")):
        pub = lab.PUBLISHED[(sistema, nome)]
        p = lab.POTENTIALS[nome](pub["p1"], pub["p2"])
        _, r0, _ = lab.scattering(p)
        d = r0 / base[(sistema, nome)][1] - 1.0
        print(f"    {eps:.0e}  {nome:6s}  {d:+.3e}")
lab.EPS_TAIL = 1e-6


# ===========================================================================
titulo(6, "As formas fechadas do PDF previem o R que o codigo usa?")
# ===========================================================================
# Eqs. 16-17 com os expoentes novos: R = (C6/(2 EPS_TAIL))^(1/4),
#                                r_min = (C12/(2 EPS_CORE))^(1/10)
# Sinal de alerta se nao baterem: a equacao do PDF nao descreve o codigo.
for sistema in ("nn", "deuteron"):
    pub = lab.PUBLISHED[(sistema, "lj")]
    p = lab.LennardJones(pub["p1"], pub["p2"])
    R_form = (pub["p1"] / (2.0 * lab.EPS_TAIL)) ** 0.25
    rm_form = (pub["p2"] / (2.0 * lab.EPS_CORE)) ** 0.1
    print(f"  {sistema:10s} R    codigo {p.R:12.6f}   formula {R_form:12.6f}"
          f"{veredito(p.R, R_form, 1e-9)}")
    print(f"  {sistema:10s} rmin codigo {p.r_min:12.6f}   formula {rm_form:12.6f}"
          f"{veredito(p.r_min, rm_form, 1e-9)}")

print("\n  e os que NAO tem forma fechada, para confirmar que brentq e mesmo")
print("  necessario (Pöschl-Teller e gaussiano com o fator r^2):")
for nome in ("mpt", "gauss"):
    pub = lab.PUBLISHED[("deuteron", nome)]
    p = lab.POTENTIALS[nome](pub["p1"], pub["p2"])
    print(f"    {nome:6s} R = {p.R:.6f}   |V(R)| R^2 = {abs(p.V(p.R))*p.R**2:.3e}"
          f"   (alvo {lab.EPS_TAIL:.0e})")


# ===========================================================================
titulo(7, "Citacoes pendentes no codigo")
# ===========================================================================
import re
fonte = open("lab.py", encoding="utf-8").read()
citadas = set(re.findall(r"\[(\d+)\]", fonte))
declaradas = {k.strip("[]") for k in lab.SOURCES if k.startswith("[")}
print(f"  citadas no arquivo:  {sorted(citadas, key=int)}")
print(f"  declaradas em SOURCES: {sorted(declaradas, key=int)}")
faltam = sorted(citadas - declaradas, key=int)
sobram = sorted(declaradas - citadas, key=int)
print(f"  CITADAS SEM FONTE: {faltam if faltam else 'nenhuma'}")
print(f"  fontes nao usadas: {sobram if sobram else 'nenhuma'}")


# ===========================================================================
titulo(8, "Quantos testes existem de fato")
# ===========================================================================
teste = open("test_lab.py", encoding="utf-8").read()
nomes = re.findall(r"^def (test_\w+)", teste, re.M)
print(f"  {len(nomes)} testes. O PDF diz 8.")
usa_solver = [n for n in nomes
              if "lab.scattering" in teste.split(f"def {n}")[1].split("\ndef ")[0]
              or "lab.bound_energy" in teste.split(f"def {n}")[1].split("\ndef ")[0]
              or "lab.tune" in teste.split(f"def {n}")[1].split("\ndef ")[0]]
print(f"\n  AUDITORIA DE CIRCULARIDADE: quais tocam o solver?")
for n in nomes:
    marca = "solver" if n in usa_solver else "SO FORMA FECHADA"
    print(f"    {marca:18s} {n}")


print("\n" + "=" * LARGURA)
print("  fim. nada foi alterado.")
print("=" * LARGURA)

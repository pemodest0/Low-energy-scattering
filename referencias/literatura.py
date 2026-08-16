# -*- coding: utf-8 -*-
"""
===============================================================================
 LITERATURA — o arquivo-mãe
===============================================================================

Todo número de literatura que este laboratório usa, num lugar só, com:

    * de qual artigo veio          (FONTES)
    * em que convenção ele vale    (CONVENCOES)
    * sob que condições é válido   (campo `valido_se` de cada registro)
    * onde a literatura discorda   (DIVERGENCIAS)

Regra da casa: **nenhum número entra aqui sem fonte**. Se a fonte não foi
conferida no PDF original, o campo `conferido` fica False e aparece marcado.

-------------------------------------------------------------------------------
COMO USAR
-------------------------------------------------------------------------------

    python referencias/literatura.py              # índice geral
    python referencias/literatura.py gauss        # tudo sobre o gaussiano
    python referencias/literatura.py deuteron     # tudo sobre o dêuteron
    python referencias/literatura.py C1a          # o que é essa convenção
    python referencias/literatura.py D1           # a divergência D1

ou, de dentro do Python:

    from referencias.literatura import ver, params, fonte, convencao, divergencia

    ver(potencial="gauss")           # tabela filtrada
    ver(sistema="dêuteron")
    params("gauss", "deuteron")      # -> (1.9102, 0.6754)
    fonte("etrych2023")              # citação completa
    convencao("C1a")                 # o que muda de convenção pra convenção
    divergencia("D1")                # onde a literatura briga consigo mesma

-------------------------------------------------------------------------------
UNIDADES — leia antes de usar qualquer número
-------------------------------------------------------------------------------

O núcleo numérico é adimensional: hbar = m_r = 1, comprimentos em unidades de
l = 1 fm. Os números daqui vêm em unidades físicas e só se comparam depois de
"recuperar as unidades" (ver src/constantes.py).

    a0 = raio de Bohr = 0.0529177 nm
    fm = 1e-15 m        angstrom = 1e-10 m       G = gauss (campo magnético)

===============================================================================
"""

from __future__ import annotations
import sys
import unicodedata


def _n(t):
    """normaliza para busca: minúsculo e sem acento (deuteron == dêuteron)."""
    t = unicodedata.normalize("NFKD", str(t))
    return "".join(c for c in t if not unicodedata.combining(c)).lower()

# =============================================================================
#  FONTES
# =============================================================================
# acesso: "aberto" = PDF livre;  "fechado" = precisa de assinatura/VPN USP
# conferido: True só se alguém abriu o PDF e leu o número na página.

FONTES = {
    # ---------------------------------------------------- o artigo reproduzido
    "rbef2023": dict(
        autores="Macêdo-Lima, M.; Madeira, L.",
        titulo="Scattering length and effective range of microscopic two-body potentials",
        revista="Rev. Bras. Ensino Fís. 45, e20230079",
        ano=2023, doi="10.1590/1806-9126-RBEF-2023-0079", acesso="aberto",
        papel="O artigo que este laboratório reproduz inteiro."),

    # ---------------------------------------------------- dois corpos, dados
    "hackenburg2006": dict(
        autores="Hackenburg, R. W.",
        titulo="Low-energy neutron-proton effective range parameters",
        revista="Phys. Rev. C 73, 044002",
        ano=2006, doi="10.1103/PhysRevC.73.044002", acesso="fechado",
        papel="Fonte de a e r0 do dêuteron. É a ref. [23] do RBEF 2023."),

    "cencek2012": dict(
        autores="Cencek, W.; Przybytek, M.; Komasa, J.; Mehl, J. B.; "
                "Jeziorski, B.; Szalewicz, K.",
        titulo="Effects of adiabatic, relativistic, and QED corrections on the "
               "pair potential of helium",
        revista="J. Chem. Phys. 136, 224303",
        ano=2012, doi="10.1063/1.4712218", acesso="fechado",
        papel="Potencial He-He ab initio. Fonte de a e r0 do dímero de 4He. "
              "Ref. [22] do RBEF 2023."),

    # ---------------------------------------------------- 39K / Feshbach
    "derrico2007": dict(
        autores="D'Errico, C.; et al.",
        titulo="Feshbach resonances in ultracold 39K",
        revista="New J. Phys. 9, 223",
        ano=2007, doi="10.1088/1367-2630/9/7/223", acesso="aberto",
        papel="Mapeamento original das ressonâncias do 39K. É o que o "
              "src/feshbach.py usa como base."),

    "zaccanti2009": dict(
        autores="Zaccanti, M.; et al.",
        titulo="Observation of an Efimov spectrum in an atomic system",
        revista="Nature Physics 5, 586",
        ano=2009, doi="10.1038/nphys1334", acesso="fechado",
        arxiv="0904.4453",
        papel="Primeiro espectro de Efimov no 39K. Ref. [14] do seu projeto."),

    "etrych2023": dict(
        autores="Etrych, J.; Martirosyan, G.; Cao, A.; Glidden, J. A. P.; "
                "Dogra, L. H.; Hutson, J. M.; Hadzibabic, Z.; Eigen, C.",
        titulo="Pinpointing Feshbach resonances and testing Efimov "
               "universalities in 39K",
        revista="Phys. Rev. Research 5, 013174",
        ano=2023, doi="10.1103/PhysRevResearch.5.013174", acesso="aberto",
        papel="Medida moderna de alta precisão no MESMO sistema do seu "
              "projeto. Discorda do Zaccanti 2009 — ver DIVERGÊNCIA D3."),

    # ---------------------------------------------------- Efimov, teoria
    "efimov1970": dict(
        autores="Efimov, V.",
        titulo="Energy levels arising from resonant two-body forces in a "
               "three-body system",
        revista="Physics Letters B 33, 563",
        ano=1970, doi="10.1016/0370-2693(70)90349-7", acesso="fechado",
        papel="O artigo original. Ref. [5] do seu projeto.",
        conferido=False),

    "braaten2006": dict(
        autores="Braaten, E.; Hammer, H.-W.",
        titulo="Universality in few-body systems with large scattering length",
        revista="Phys. Rep. 428, 259",
        ano=2006, doi="10.1016/j.physrep.2006.03.001", acesso="fechado",
        papel="A revisão canônica. Fonte de s0 e das razões universais. "
              "Ref. [6] do seu projeto."),

    "naidon2017": dict(
        autores="Naidon, P.; Endo, S.",
        titulo="Efimov physics: a review",
        revista="Rep. Prog. Phys. 80, 056001",
        ano=2017, doi="10.1088/1361-6633/aa50e8", acesso="fechado",
        papel="Revisão moderna; de onde vem a Fig. 1 do seu projeto. "
              "Ref. [7] do seu projeto.",
        conferido=False),

    "chin2010": dict(
        autores="Chin, C.; Grimm, R.; Julienne, P.; Tiesinga, E.",
        titulo="Feshbach resonances in ultracold gases",
        revista="Rev. Mod. Phys. 82, 1225",
        ano=2010, doi="10.1103/RevModPhys.82.1225", acesso="fechado",
        papel="A referência para tudo que é Feshbach: a_bg, Delta, s_res, R_vdW."),

    # ---------------------------------------------------- o grupo
    "madeira2021": dict(
        autores="Madeira, L.; Frederico, T.; Gandolfi, S.; Tomio, L.; "
                "Yamashita, M. T.",
        titulo="Quantum Monte Carlo studies of a trimer scaling function with "
               "microscopic two- and three-body interactions",
        revista="Phys. Rev. A 104, 033301",
        ano=2021, doi="10.1103/PhysRevA.104.033301", acesso="fechado",
        papel="*** O artigo que a sua dissertação estende para sistemas "
              "confinados. Ref. [18] do seu projeto. ***"),

    "madeira2024": dict(
        autores="Madeira, L.",
        titulo="The Role of the Effective Range in Strongly-Interacting "
               "Few-Body Systems",
        revista="Few-Body Systems 65, 70",
        ano=2024, doi="10.1007/s00601-024-01940-2", acesso="fechado",
        papel="*** A ponte entre o r0 de dois corpos e a física de três "
              "corpos. É por isso que o alcance efetivo veio primeiro. "
              "Ref. [29] do seu projeto. ***",
        conferido=False),

    # ---------------------------------------------------- benchmarks externos
    "aziz_hfdb": dict(
        autores="Aziz, R. A.; et al.",
        titulo="Potencial He-He HFD-B",
        revista="PENDENTE — é a ref. [25] do RBEF 2023",
        ano=None, doi="PENDENTE", acesso="?",
        papel="Potencial realista de hélio. Liga o dímero; o Lennard-Jones "
              "clássico não liga. Ver ACHADO 3.",
        conferido=False),

    "jeszenszki": dict(
        autores="Jeszenszki, P.; et al.",
        titulo="a(v) do potencial gaussiano em forma fechada",
        revista="PENDENTE — é a ref. [34] do RBEF 2023",
        ano=None, doi="PENDENTE", acesso="?",
        papel="Fórmula fechada para o gaussiano. Nosso benchmark externo mais "
              "preciso (concorda em ~1e-6).",
        conferido=False),

    "bruch_tjon": dict(
        autores="Bruch, L. W.; Tjon, J. A.",
        titulo="Three-body system in two dimensions",
        revista="PENDENTE — conferir",
        ano=None, doi="PENDENTE", acesso="?",
        papel="Em 2D há exatamente DOIS trímeros universais, não uma torre. "
              "É a âncora do lado 2D do seu crossover dimensional.",
        conferido=False),
}


# =============================================================================
#  CONVENÇÕES — a pedra de Roseta
# =============================================================================
# Dois artigos podem publicar números diferentes para a MESMA física só porque
# escreveram o prefator do potencial de outro jeito. Antes de comparar
# qualquer número com qualquer outro, veja se a convenção bate.

CONVENCOES = {
    "C1a": dict(
        titulo="Potencial com hbar^2/(2 m_r)",
        formula="V(r) = (hbar^2 / (2 m_r)) * U(r)   =>   u'' = U u",
        nota="É a convenção em que a Tabela 4 (Lennard-Jones) do RBEF fecha. "
             "Ver DIVERGÊNCIA D1."),

    "C1b": dict(
        titulo="Potencial com hbar^2/m_r  (a do nosso código)",
        formula="V(r) = (hbar^2 / m_r) * U(r)       =>   u'' = 2 V u",
        nota="É a convenção do src/. As Tabelas 2 e 3 do RBEF fecham aqui."),

    "C2a": dict(
        titulo="Escala de van der Waals",
        formula="R_vdW = (1/2) (2 m_r C6 / hbar^2)^(1/4)",
        nota="Cuidado: alguns artigos chamam de l_vdW o dobro disso. "
             "Para o 39K, R_vdW = 64.6 a0 (etrych2023)."),

    "C3": dict(
        titulo="Sinal do comprimento de espalhamento",
        formula="u(r) ~ r - a  fora do alcance;  a > 0 <=> existe estado ligado raso",
        nota="Convenção de física atômica/nuclear moderna. Textos antigos às "
             "vezes usam o sinal oposto."),

    "C4a": dict(
        titulo="Expansão de alcance efetivo",
        formula="k cot(delta_0) = -1/a + (1/2) r0 k^2 + O(k^4)",
        nota="O fator 1/2 é o padrão. Alguns textos absorvem o 1/2 no r0 e "
             "publicam metade do valor."),

    "C5": dict(
        titulo="Unidades de energia",
        formula="MeV para física nuclear; mK ou nK para átomos frios; G para campo",
        nota="1 mK = 1e-3 K. E = -hbar^2/(2 m_r a^2) só vale no limite de "
             "alcance zero."),

    "C6": dict(
        titulo="Parâmetro de três corpos",
        formula="kappa_* (número de onda na unitariedade) ou a_- (onde o "
                "trímero cruza o limiar)",
        nota="São o MESMO conteúdo físico em roupas diferentes. Artigos "
             "experimentais publicam a_-; teóricos publicam kappa_*."),

    "C7": dict(
        titulo="Razões de Efimov",
        formula="s0 = 1.0062378;  E_n/E_(n+1) = e^(2pi/s0) = 515.03;  "
                "razão de escala = e^(pi/s0) = 22.694",
        nota="Só valem para TRÊS BÓSONS IDÊNTICOS em alcance zero e 3D. "
             "Massas diferentes mudam s0; dimensão diferente também."),

    "C8": dict(
        titulo="Alcance do potencial vs alcance efetivo",
        formula="R = onde V(R) ~ 0   !=   r0 = 2 int [g0^2 - u0^2] dr",
        nota="São coisas diferentes e às vezes têm valores parecidos. "
             "Não confunda: R é geometria do potencial, r0 é da solução."),
}


# =============================================================================
#  POTENCIAIS — a forma funcional, e quando cada um vale a pena
# =============================================================================

POTENCIAIS = {
    "poco": dict(
        nome="Poço esférico",
        formula="V(r) = -v mu^2  (r < R);  0  (r >= R);   R = 1/mu",
        parametros=("v", "mu"),
        fonte="rbef2023", local="Eq. (74)",
        a_analitico="a = R [1 - tan(x)/x],  x = sqrt(2v)",
        limiar="v > pi^2/8 = 1.2337  para o 1o estado ligado",
        bom_para="Ter resposta analítica exata. É a âncora de validação nº 1.",
        cuidado="DESCONTÍNUO na borda. Numerov perde a ordem alta aqui e "
                "degrada de O(dr^4) para O(dr^1) — ver ACHADO 2."),

    "mpt": dict(
        nome="Pöschl-Teller modificado",
        formula="V(r) = -v mu^2 / cosh^2(mu r)",
        parametros=("v", "mu"),
        fonte="rbef2023", local="Eq. (116)-(117)",
        a_analitico="a mu = (pi/2) cot(pi lam/2) + gamma + Psi(lam),  "
                    "lam = (1 + sqrt(1+8v))/2",
        limiar="v = 1 é a unitariedade exata (lam = 2), com r0 = 2/mu",
        bom_para="Ser LISO e ter resposta analítica. O melhor caso de teste "
                 "para medir a ordem de um integrador.",
        cuidado="Cauda exponencial: o 'alcance' depende do corte que você "
                "escolher (EPS_CORTE no código)."),

    "gauss": dict(
        nome="Gaussiano",
        formula="V(r) = -v mu^2 exp(-r^2 mu^2)",
        parametros=("v", "mu"),
        fonte="rbef2023", local="Eq. (120)",
        a_analitico="forma fechada em jeszenszki (ref. [34] do RBEF)",
        limiar="v ~ 1.342 para mu = 1",
        bom_para="Liso, decai rápido, e tem benchmark externo independente "
                 "(concordância ~1e-6).",
        cuidado=None),

    "lj": dict(
        nome="Lennard-Jones",
        formula="V(r) = (1/2) (C12/r^12 - C6/r^6)     <-- note o 1/2",
        parametros=("C12", "C6"),
        fonte="rbef2023", local="Eq. (121) + Tabela 4",
        a_analitico=None,
        limiar=None,
        bom_para="Ser realista: caroço repulsivo + cauda de van der Waals. "
                 "É a forma que descreve átomos de verdade.",
        cuidado="*** O fator 1/2 NÃO está na Eq. (121) publicada. Sem ele a "
                "Tabela 4 não fecha. Ver DIVERGÊNCIA D1. *** "
                "Além disso, diverge em r->0: precisa de caroço numérico "
                "(r_min tal que U(r_min) ~ 1e10)."),

    "aziz_hfdb": dict(
        nome="Aziz HFD-B (He-He)",
        formula="forma HFD com termos de dispersão C6, C8, C10",
        parametros=("tabelados",),
        fonte="aziz_hfdb", local="ref. [25] do RBEF",
        a_analitico=None, limiar=None,
        bom_para="Hélio de verdade. Liga o dímero (a = +88.4 A).",
        cuidado="*** O Lennard-Jones clássico de de Boer NÃO liga o dímero "
                "(a = -178 A). Mesmo átomo, respostas opostas. "
                "Ver ACHADO 3. ***"),
}


# =============================================================================
#  VALORES — o corpo do arquivo
# =============================================================================
# Cada registro:
#   grandeza, sistema, potencial, params, valor, unidade, tipo, fonte,
#   local (onde no artigo), convencao, valido_se, nota, conferido

def _v(**kw):
    kw.setdefault("params", None); kw.setdefault("potencial", None)
    kw.setdefault("nota", None); kw.setdefault("conferido", True)
    kw.setdefault("incerteza", None); kw.setdefault("local", None)
    return kw

VALORES = [

    # ---------- Tabela 3 do RBEF: parâmetros ajustados ---------------------
    _v(grandeza="a", sistema="nêutron-nêutron", potencial="poco",
       params=(1.1096, 0.3918), valor=-18.52, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C1b",
       valido_se="unidades adimensionais do artigo; onda-s; energia zero"),
    _v(grandeza="r0", sistema="nêutron-nêutron", potencial="poco",
       params=(1.1096, 0.3918), valor=2.70, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C4a",
       valido_se="idem"),

    _v(grandeza="a", sistema="nêutron-nêutron", potencial="mpt",
       params=(0.9071, 0.7991), valor=-18.51, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C1b",
       valido_se="idem"),
    _v(grandeza="a", sistema="nêutron-nêutron", potencial="gauss",
       params=(1.2121, 0.5672), valor=-18.55, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C1b",
       valido_se="idem"),

    _v(grandeza="a", sistema="dêuteron (n-p tripleto)", potencial="poco",
       params=(1.7575, 0.5000), valor=5.40, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C1b",
       valido_se="idem"),
    _v(grandeza="a", sistema="dêuteron (n-p tripleto)", potencial="mpt",
       params=(1.4388, 0.8631), valor=5.40, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C1b",
       valido_se="idem", nota="r0 publicado = 1.73 (os outros dão 1.70)"),
    _v(grandeza="a", sistema="dêuteron (n-p tripleto)", potencial="gauss",
       params=(1.9102, 0.6754), valor=5.40, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C1b",
       valido_se="idem"),
    _v(grandeza="r0", sistema="dêuteron (n-p tripleto)", potencial="gauss",
       params=(1.9102, 0.6754), valor=1.70, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C4a",
       valido_se="idem"),

    _v(grandeza="r0", sistema="unitariedade (|a| -> inf)", potencial="gauss",
       params=(1.3420, 1.4349), valor=1.00, unidade="fm", tipo="ajustado",
       fonte="rbef2023", local="Tabela 3", convencao="C4a",
       valido_se="1/a = 0 exatamente"),
    _v(grandeza="v_limiar", sistema="unitariedade", potencial="poco",
       params=(1.2337, 1.0), valor=1.2337005501, unidade="adimensional",
       tipo="teorico", fonte="rbef2023", local="Eq. (86)", convencao="C1b",
       valido_se="pi^2/8 exato", nota="onde nasce o 1o estado ligado do poço"),

    # ---------- Tabela 4 do RBEF: Lennard-Jones (ATENÇÃO à convenção) ------
    _v(grandeza="a", sistema="nêutron-nêutron", potencial="lj",
       params=(3.08836698, 9.86668911), valor=-18.5, unidade="fm",
       tipo="ajustado", fonte="rbef2023", local="Tabela 4", convencao="C1a",
       valido_se="SÓ em C1a. A Eq. (121) impressa está em C1b e não fecha.",
       nota="DIVERGÊNCIA D1 — params na ordem (C12, C6)"),
    _v(grandeza="r0", sistema="nêutron-nêutron", potencial="lj",
       params=(3.08836698, 9.86668911), valor=2.71, unidade="fm",
       tipo="ajustado", fonte="rbef2023", local="Tabela 4", convencao="C1a",
       valido_se="idem — ver D1"),
    _v(grandeza="a", sistema="dêuteron (n-p tripleto)", potencial="lj",
       params=(0.90485319, 6.81472000), valor=5.4, unidade="fm",
       tipo="ajustado", fonte="rbef2023", local="Tabela 4", convencao="C1a",
       valido_se="idem — ver D1"),

    # ---------- Tabela 1 do RBEF: valores experimentais/ab initio ----------
    _v(grandeza="a", sistema="dêuteron (n-p tripleto, 3S1)", valor=5.4112,
       unidade="fm", tipo="experimental", fonte="hackenburg2006",
       local="citado na Tabela 1 do RBEF", convencao="C3",
       valido_se="canal tripleto 3S1; baixa energia",
       nota="PDF fechado — valor e incerteza ainda não conferidos na fonte",
       conferido=False),
    _v(grandeza="r0", sistema="dêuteron (n-p tripleto, 3S1)", valor=1.7436,
       unidade="fm", tipo="experimental", fonte="hackenburg2006",
       convencao="C4a", valido_se="idem", conferido=False),
    _v(grandeza="E_b", sistema="dêuteron", valor=-2.224, unidade="MeV",
       tipo="experimental", fonte="rbef2023", local="Tabela 1", convencao="C5",
       valido_se="energia de ligação medida"),
    _v(grandeza="E_zr", sistema="dêuteron", valor=-1.416, unidade="MeV",
       tipo="teorico", fonte="rbef2023", local="Tabela 1 / Eq. (93)",
       convencao="C5", valido_se="alcance zero: E = -hbar^2/(2 m_r a^2)",
       nota="erra 36% — mostra que r0 NÃO é detalhe no dêuteron"),
    _v(grandeza="E_fr", sistema="dêuteron", valor=-2.223, unidade="MeV",
       tipo="teorico", fonte="rbef2023", local="Tabela 1 / Eq. (95)",
       convencao="C5", valido_se="alcance finito: 1/a = kappa - r0 kappa^2/2",
       nota="acerta em 0.05% — o r0 vale 0.8 MeV aqui"),

    _v(grandeza="a", sistema="dímero de 4He", valor=90.4, unidade="angstrom",
       tipo="teorico", fonte="cencek2012", local="citado na Tabela 1 do RBEF",
       convencao="C3", valido_se="potencial ab initio com correções QED",
       nota="PDF fechado — não conferido na fonte", conferido=False),
    _v(grandeza="r0", sistema="dímero de 4He", valor=8.0, unidade="angstrom",
       tipo="teorico", fonte="cencek2012", convencao="C4a",
       valido_se="idem", conferido=False),
    _v(grandeza="E_b", sistema="dímero de 4He", valor=-1.62, unidade="mK",
       tipo="teorico", fonte="cencek2012", convencao="C5",
       valido_se="idem", conferido=False),

    _v(grandeza="a", sistema="dímero de 4He", potencial="aziz_hfdb",
       valor=88.43, unidade="angstrom", tipo="calculado (nosso)",
       fonte="aziz_hfdb", convencao="C3",
       valido_se="Aziz HFD-B; nosso integrador",
       nota="LIGA o dímero (1 nó). Compare com o LJ logo abaixo."),
    _v(grandeza="a", sistema="dímero de 4He", potencial="lj",
       valor=-177.95, unidade="angstrom", tipo="calculado (nosso)",
       fonte="rbef2023", convencao="C3",
       valido_se="Lennard-Jones clássico de de Boer",
       nota="*** ACHADO 3: NÃO liga o dímero (0 nós). Mesmo átomo, sinal "
            "oposto de a. A forma do potencial decide. ***"),

    # ---------- benchmarks externos ----------------------------------------
    _v(grandeza="a", sistema="par n-p singleto", valor=-23.74, unidade="fm",
       tipo="experimental", fonte="hackenburg2006", convencao="C3",
       valido_se="canal singleto 1S0",
       nota="estado VIRTUAL, não ligado (kappa < 0). |a| enorme com potencial "
            "curto — o exemplo de livro de 'quase liga'.", conferido=False),
    _v(grandeza="r0", sistema="par n-p singleto", valor=2.77, unidade="fm",
       tipo="experimental", fonte="hackenburg2006", convencao="C4a",
       valido_se="idem", conferido=False),

    _v(grandeza="r0", sistema="esfera dura", valor=None, unidade="—",
       tipo="teorico", fonte="rbef2023", convencao="C4a",
       valido_se="limite de esfera dura de raio R",
       nota="r0 = 2R/3 exato. Nosso teste de esfera mole converge para isso."),

    # ---------- Efimov: constantes universais ------------------------------
    _v(grandeza="s0", sistema="3 bósons idênticos, alcance zero, 3D",
       valor=1.0062378, unidade="adimensional", tipo="teorico",
       fonte="braaten2006", convencao="C7",
       valido_se="TRÊS BÓSONS IDÊNTICOS, alcance zero, unitariedade, d=3. "
                 "Massas diferentes ou dimensão diferente mudam s0.",
       nota="raiz de  s cosh(pi s/2) = (8/sqrt3) sinh(pi s/6)"),
    _v(grandeza="E_n/E_(n+1)", sistema="3 bósons idênticos, unitariedade",
       valor=515.03, unidade="adimensional", tipo="teorico",
       fonte="braaten2006", convencao="C7", valido_se="idem",
       nota="= e^(2 pi / s0)"),
    _v(grandeza="razão de escala", sistema="3 bósons idênticos, unitariedade",
       valor=22.694, unidade="adimensional", tipo="teorico",
       fonte="braaten2006", convencao="C7", valido_se="idem",
       nota="= e^(pi / s0)"),

    _v(grandeza="E3/E2 (fundamental)", sistema="3 bósons em 2D",
       valor=16.522, unidade="adimensional", tipo="teorico",
       fonte="bruch_tjon", convencao="C7",
       valido_se="DUAS dimensões, alcance zero",
       nota="*** Em 2D há exatamente DOIS trímeros, não uma torre. Esta é a "
            "sua âncora do lado 2D do crossover. FONTE A CONFERIR. ***",
       conferido=False),
    _v(grandeza="E3/E2 (excitado)", sistema="3 bósons em 2D",
       valor=1.2704, unidade="adimensional", tipo="teorico",
       fonte="bruch_tjon", convencao="C7", valido_se="idem",
       conferido=False),

    # ---------- 39K: Feshbach ----------------------------------------------
    _v(grandeza="B_res", sistema="39K |1,1>, ressonância larga", valor=402.74,
       incerteza=0.01, unidade="G", tipo="experimental", fonte="etrych2023",
       local="Tabela I", convencao="C5",
       valido_se="estado hiperfino |1,1>; medida de alta precisão 2023",
       nota="Nosso src/feshbach.py ainda usa 402.50 (Zaccanti 2009). "
            "DECISÃO PENDENTE: atualizar?"),
    _v(grandeza="B_res", sistema="39K |1,1>, ressonância larga", valor=402.4,
       unidade="G", tipo="teorico", fonte="derrico2007",
       local="Tabela 1, coluna B_th", convencao="C5",
       valido_se="cálculo de canais acoplados de 2007",
       nota="B_exp da época = 403.4(7) G. Superado por etrych2023."),
    _v(grandeza="B_zero", sistema="39K |1,1>, cruzamento de zero", valor=350.4,
       incerteza=0.1, unidade="G", tipo="experimental", fonte="etrych2023",
       local="Tabela I", convencao="C5", valido_se="idem"),
    _v(grandeza="a_bg*Delta", sistema="39K |1,1>", valor=1530, incerteza=20,
       unidade="a0*G", tipo="experimental", fonte="etrych2023",
       local="Tabela I", convencao="C5", valido_se="idem",
       nota="nosso código dá 1511 a0.G — compatível"),
    _v(grandeza="s_res", sistema="39K |1,1>, 402.74 G", valor=2.8,
       unidade="adimensional", tipo="teorico", fonte="etrych2023",
       local="Tabela III", convencao="C5",
       valido_se="força da ressonância",
       nota="INTERMEDIÁRIA — nem larga nem estreita. É a chave da D2."),
    _v(grandeza="r_eff", sistema="39K |1,1>, perto de 402.74 G", valor=136,
       unidade="a0", tipo="teorico", fonte="etrych2023",
       local="Tabela III (canais acoplados)", convencao="C4a",
       valido_se="perto da ressonância"),
    _v(grandeza="R_vdW", sistema="39K", valor=64.6, unidade="a0",
       tipo="teorico", fonte="etrych2023", local="Eq. (2)", convencao="C2a",
       valido_se="escala de van der Waals do K-K",
       nota="confirmado independentemente — nosso valor está certo"),
    _v(grandeza="a_bar", sistema="39K", valor=61.76, unidade="a0",
       tipo="teorico", fonte="etrych2023", convencao="C2a",
       valido_se="comprimento médio de espalhamento"),

    # ---------- 39K: Efimov ------------------------------------------------
    _v(grandeza="a_menos", sistema="39K, universalidade vdW", valor=-630,
       unidade="a0", tipo="teorico", fonte="etrych2023", convencao="C6",
       valido_se="universalidade de Efimov-van der Waals: a_- ~ -9.7 R_vdW",
       nota="ver D2 (RESOLVIDA)"),
    _v(grandeza="a_menos", sistema="39K, ressonância de 402.74 G", valor=-830,
       incerteza=40, unidade="a0", tipo="experimental", fonte="etrych2023",
       convencao="C6", valido_se="medida na ressonância intermediária",
       nota="DIVERGE do Zaccanti 2009 (-1500) na MESMA ressonância — ver D3"),
    _v(grandeza="a_menos", sistema="39K, ressonância de ~402 G", valor=-1500,
       unidade="a0", tipo="experimental", fonte="zaccanti2009",
       convencao="C6", valido_se="medida de 2009",
       nota="*** DIVERGÊNCIA D3: Etrych 2023 mede -830(40) no mesmo sistema. "
            "Fator ~1.8 de diferença. ***"),
    _v(grandeza="eta_menos", sistema="39K", valor=0.14, unidade="adimensional",
       tipo="experimental", fonte="zaccanti2009", convencao="C6",
       valido_se="largura inelástica",
       nota="Etrych mede 0.27(6) — parte da D3"),
    _v(grandeza="a_p/a_menos", sistema="39K, todas as ressonâncias",
       valor=-1.08, unidade="adimensional", tipo="experimental",
       fonte="etrych2023", convencao="C6", valido_se="teste de universalidade"),
]


# =============================================================================
#  DIVERGÊNCIAS — onde a literatura briga consigo mesma
# =============================================================================
# Regra: divergência NÃO se apaga. Registra-se, com as duas versões.

DIVERGENCIAS = {
    "D1": dict(
        titulo="Fator 2 na Eq. (121) do RBEF 2023 (Lennard-Jones)",
        status="RESOLVIDA (nosso lado)",
        o_que="A Eq. (121) impressa escreve V = (hbar^2/m_r)[C12/r^12 - C6/r^6], "
              "mas as constantes da Tabela 4 SÓ reproduzem (a, r0) se o "
              "prefator for hbar^2/(2 m_r).",
        como_achamos="Reproduzindo a Tabela 4 e vendo que não fechava. Com o "
                     "fator 1/2 fecha em ~0.1%.",
        acao="src/potenciais.py usa V = (1/2)(C12/r^12 - C6/r^6). Documentado "
             "na classe LennardJones.",
        fontes=["rbef2023"]),

    "D2": dict(
        titulo="a_menos do 39K vs universalidade de van der Waals",
        status="RESOLVIDA",
        o_que="A universalidade Efimov-vdW prevê a_- ~ -9.7 R_vdW = -630 a0, "
              "mas o experimento mede -830(40) a0.",
        como_achamos="Comparando etrych2023 com a previsão universal.",
        acao="A ressonância de 402.74 G tem s_res = 2.8 — força INTERMEDIÁRIA. "
             "A universalidade vdW pressupõe ressonância larga (s_res >> 1). "
             "Não é falha da teoria, é hipótese fora de validade.",
        fontes=["etrych2023"]),

    "D3": dict(
        titulo="Zaccanti 2009 vs Etrych 2023 — dois experimentos, uma ressonância",
        status="ABERTA",
        o_que="Zaccanti mede a_- = -1500 a0 e eta_- = 0.14. Etrych mede "
              "a_- = -830(40) a0 e eta_- = 0.27(6). Mesmo átomo, mesma "
              "ressonância, fator ~1.8 de diferença.",
        como_achamos="Ao montar o BENCHMARKS.yaml e colocar os dois lado a lado.",
        acao="NÃO escolher um. Registrar os dois e usar o mais recente com "
             "ressalva. Se a sua dissertação encostar em a_- do 39K, isto "
             "precisa ser discutido — e é assunto para o Lucas e a Patrícia.",
        fontes=["zaccanti2009", "etrych2023"]),
}


# =============================================================================
#  ACHADOS — o que este laboratório descobriu sozinho
# =============================================================================

ACHADOS = {
    1: dict(titulo="Fator 2 na Eq. (121)",
            resumo="Ver DIVERGÊNCIA D1.",
            onde="src/potenciais.py :: LennardJones"),
    2: dict(titulo="Numerov perde a ordem em potencial descontínuo",
            resumo="No poço esférico (borda abrupta) a diferença central é "
                   "O(dr^2) e o Numerov degrada para O(dr^1) — chegando a "
                   "errar ~1000x MAIS que o método simples em dr=1e-3. "
                   "A ordem alta do Numerov pressupõe V suave (classe C4).",
            onde="tests/ + notebooks/LABORATORIO_DO_ZERO.ipynb Ex 14"),
    3: dict(titulo="O dímero de hélio vive no fio da navalha",
            resumo="Lennard-Jones de de Boer: a = -178 A, NÃO liga. "
                   "Aziz HFD-B: a = +88.4 A, E = -1.69 mK, LIGA. "
                   "Mesmo átomo, sinal oposto de a. Quando |a| >> alcance, "
                   "a resposta é hipersensível à forma do potencial.",
            onde="referencias/bench_aziz.csv e bench_he_dimer.csv"),
}


# =============================================================================
#  CONSULTA
# =============================================================================

_L = 78

def _hr(c="-"): return c * _L

def fonte(fid: str):
    """Citação completa de uma fonte."""
    f = FONTES.get(fid)
    if not f:
        print(f"fonte '{fid}' não existe. Disponíveis: {', '.join(sorted(FONTES))}")
        return
    print(_hr("=")); print(f" {fid}"); print(_hr("="))
    print(f"  {f['autores']}")
    print(f"  {f['titulo']}")
    print(f"  {f['revista']} ({f['ano']})" if f["ano"] else f"  {f['revista']}")
    print(f"  doi: {f['doi']}    acesso: {f['acesso']}")
    if f.get("arxiv"): print(f"  arXiv: {f['arxiv']}")
    if not f.get("conferido", True):
        print("  [!] NÃO CONFERIDO na fonte original")
    print(f"\n  papel: {f['papel']}")
    print(_hr("=")); print()


def convencao(cid: str):
    """O que significa uma convenção, e o que muda se você trocar."""
    c = CONVENCOES.get(cid)
    if not c:
        print(f"convenção '{cid}' não existe. Disponíveis: {', '.join(sorted(CONVENCOES))}")
        return
    print(_hr("=")); print(f" {cid} · {c['titulo']}"); print(_hr("="))
    print(f"  {c['formula']}")
    print(f"\n  {c['nota']}")
    print(_hr("=")); print()


def divergencia(did: str):
    """Onde a literatura discorda de si mesma."""
    d = DIVERGENCIAS.get(did)
    if not d:
        print(f"divergência '{did}' não existe. Disponíveis: {', '.join(sorted(DIVERGENCIAS))}")
        return
    print(_hr("=")); print(f" {did} · {d['titulo']}   [{d['status']}]"); print(_hr("="))
    for rot, ch in (("o que é", "o_que"), ("como achamos", "como_achamos"),
                    ("o que fizemos", "acao")):
        print(f"\n  {rot.upper()}:")
        for ln in _quebra(d[ch], _L - 6): print(f"    {ln}")
    print(f"\n  fontes: {', '.join(d['fontes'])}")
    print(_hr("=")); print()


def potencial(pid: str):
    """Ficha completa de um potencial."""
    p = POTENCIAIS.get(pid)
    if not p:
        print(f"potencial '{pid}' não existe. Disponíveis: {', '.join(sorted(POTENCIAIS))}")
        return
    print(_hr("=")); print(f" {pid} · {p['nome']}"); print(_hr("="))
    print(f"  {p['formula']}")
    print(f"  parâmetros: {p['parametros']}")
    print(f"  fonte: {p['fonte']}" + (f" ({p['local']})" if p.get("local") else ""))
    if p.get("a_analitico"): print(f"\n  a analítico: {p['a_analitico']}")
    if p.get("limiar"):      print(f"  limiar:      {p['limiar']}")
    print(f"\n  BOM PARA: ", end="")
    print("\n            ".join(_quebra(p["bom_para"], _L - 14)))
    if p.get("cuidado"):
        print(f"\n  CUIDADO:  ", end="")
        print("\n            ".join(_quebra(p["cuidado"], _L - 14)))
    print(_hr("=")); print()


def params(pot: str, sistema_chave: str):
    """Devolve os parâmetros publicados. Ex: params('gauss', 'deuteron')."""
    alvo = _n(sistema_chave)
    for r in VALORES:
        if r["potencial"] == pot and r["params"] and alvo in _n(r["sistema"]):
            return r["params"]
    return None


def ver(potencial=None, sistema=None, grandeza=None, fonte=None):
    """Tabela filtrada. Sem argumentos, mostra tudo."""
    sel = [r for r in VALORES
           if (potencial is None or r["potencial"] == potencial)
           and (sistema is None or _n(sistema) in _n(r["sistema"]))
           and (grandeza is None or r["grandeza"] == grandeza)
           and (fonte is None or r["fonte"] == fonte)]
    if not sel:
        print("nada encontrado com esse filtro."); return
    print(_hr("="))
    print(f" {'grandeza':<18}{'valor':>13} {'un':<10}{'pot':<7}{'fonte':<15}conv")
    print(_hr())
    sist = None
    for r in sel:
        if r["sistema"] != sist:
            sist = r["sistema"]; print(f"\n  # {sist}")
        val = "—" if r["valor"] is None else f"{r['valor']:g}"
        if r["incerteza"]: val += f"({r['incerteza']:g})"
        mark = "" if r["conferido"] else "  [!]"
        print(f" {r['grandeza']:<18}{val:>13} {r['unidade']:<10}"
              f"{(r['potencial'] or '—'):<7}{r['fonte']:<15}{r['convencao']}{mark}")
        if r.get("nota"):
            for ln in _quebra(r["nota"], _L - 8): print(f"        > {ln}")
    print(_hr("="))
    n_nc = sum(1 for r in sel if not r["conferido"])
    if n_nc: print(f" [!] {n_nc} valor(es) ainda não conferido(s) na fonte original")
    print(f" {len(sel)} registro(s).  Detalhes: fonte('id'), convencao('Cx'), "
          f"potencial('nome')")
    print()


def _quebra(t, w):
    out, linha = [], ""
    for p in str(t).split():
        if len(linha) + len(p) + 1 > w:
            out.append(linha); linha = p
        else:
            linha = (linha + " " + p).strip()
    if linha: out.append(linha)
    return out or [""]


def indice():
    print(_hr("="))
    print(" LITERATURA — arquivo-mãe do laboratório de espalhamento")
    print(_hr("="))
    print(f"\n  {len(VALORES)} valores · {len(FONTES)} fontes · "
          f"{len(POTENCIAIS)} potenciais · {len(CONVENCOES)} convenções · "
          f"{len(DIVERGENCIAS)} divergências")

    print("\n  POTENCIAIS")
    for k, p in POTENCIAIS.items():
        flag = "  <-- cuidado" if p.get("cuidado") else ""
        print(f"    {k:<12} {p['nome']}{flag}")

    print("\n  SISTEMAS")
    for s in dict.fromkeys(r["sistema"] for r in VALORES):
        print(f"    {s}")

    print("\n  DIVERGÊNCIAS")
    for k, d in DIVERGENCIAS.items():
        print(f"    {k}  [{d['status']:<16}] {d['titulo']}")

    print("\n  ACHADOS DESTE LABORATÓRIO")
    for k, a in ACHADOS.items():
        print(f"    {k}. {a['titulo']}")

    nc = [r for r in VALORES if not r["conferido"]]
    if nc:
        print(f"\n  [!] {len(nc)} valores ainda NÃO conferidos no PDF original:")
        for r in nc:
            print(f"      {r['grandeza']:<20} {r['sistema']:<32} ({r['fonte']})")

    print("\n" + _hr())
    print("  uso:  python referencias/literatura.py <potencial|sistema|Cx|Dx>")
    print(_hr("=")); print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        indice()
    else:
        q = sys.argv[1]
        if q in CONVENCOES:      convencao(q)
        elif q in DIVERGENCIAS:  divergencia(q)
        elif q in POTENCIAIS:    potencial(q); ver(potencial=q)
        elif q in FONTES:        fonte(q); ver(fonte=q)
        else:                    ver(sistema=q)

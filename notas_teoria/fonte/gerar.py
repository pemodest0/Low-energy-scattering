"""Every generated number in one place: the tables, the figures and the data.

    python gerar.py              # about 20 min
    CACHE=0 python gerar.py      # ignore what is on disk and redo it all

It writes three things, all from one set of computations:

    tabelas/gerado.tex    the table bodies and every loose number the notes
                          and the write-up quote, as LaTeX macros
    ../../resultados/     the csv files and the three figures
    stdout                the same numbers, to read while it runs

This used to be two scripts -- one for the notes, one for the deliverable --
doing the same calculations twice on two different meshes. One mesh, one run.

Only the ROWS of each table are emitted, under a macro name, never the header:
the notes are in English and the write-up in Portuguese, and both take the same
rows from here and supply their own column labels. Shape names go out as
\\FormaWell and friends for the same reason.
"""
import csv
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "..", "notebooks", "2_dois_corpos"))

import etrych
import lab

SAIDA = os.path.join(AQUI, "..", "..", "resultados")
CACHE = os.path.join(AQUI, ".cache_gerar")
USAR_CACHE = os.environ.get("CACHE", "1") != "0"

FORMAS = ("well", "mpt", "gauss", "lj")
MACRO = {"well": r"\FormaWell", "mpt": r"\FormaMpt",
         "gauss": r"\FormaGauss", "lj": r"\FormaLj"}
COR = {"well": "#2471a3", "mpt": "#c0392b", "gauss": "#1e8449", "lj": "#8e44ad"}

# One mesh, 31 points. The step of 0.02 lands exactly on -0.16, the far end of
# the band where the five r_eff/a_- sit, and on 0.0, the unitary seed. The
# endpoints reach past the positive feature at r_eff/a_p ~ +0.16. Wider was
# measured -- out to -0.90 and +0.45, where the shapes stop, the positive end
# because the effective-range expansion has no bound state past r0/a = 1/2 --
# and the curves stay straight out there, so the extra range is picture and
# not result.
ALVOS = [round(-0.30 + i * 0.02, 10) for i in range(31)]

SAIDA_TEX = []          # (macro name, lines)


def memo(nome, funcao):
    """Run `funcao`, or read back what it produced the last time it ran.

    The run is long and made of independent pieces, so an interrupted one
    resumes instead of starting over. The cache holds only what this script
    computes from lab.py and etrych.py, so deleting it costs time and nothing
    else.
    """
    caminho = os.path.join(CACHE, nome + ".json")
    if USAR_CACHE and os.path.exists(caminho):
        print("  reused", nome)
        return json.load(open(caminho, encoding="utf-8"))
    valor = funcao()
    if not os.path.isdir(CACHE):
        os.makedirs(CACHE)
    json.dump(valor, open(caminho, "w", encoding="utf-8"))
    return valor


def notacao(x):
    """A number for LaTeX: plain if it reads well, scientific if it does not."""
    if x == 0.0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    if -4 < e < 6:
        return "%g" % x
    return r"\ensuremath{%g \times 10^{%d}}" % (x / 10.0 ** e, e)


# =============================================================================
#  What gets computed
# =============================================================================
def calibracoes():
    """Both calibrations, five rows, four shapes. Forty records."""
    fora = []
    for row in etrych.TABLE_III:
        for name in FORMAS:
            for rotulo, par in (("unitaria", etrych.unitary(name, row["r_eff"])),
                                ("a_menos", etrych.calibrate(name, row))):
                p1, p2 = par
                a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
                fora.append({"calibracao": rotulo, "estado": row["state"],
                             "B_res": row["B_txt"], "s_res": row["s_res"],
                             "r_eff": row["r_eff"], "a_menos": row["a_minus"],
                             "forma": name, "p1": p1, "p2": p2,
                             "a": a, "inv_a": 1.0 / a, "r0": r0,
                             "estados": nodes})
    return fora


def familias():
    """The two prescriptions for every shape, plus the C6-alone walk."""
    fora = {}
    for name in FORMAS:
        fora[name + "|reajustada"] = memo(
            "fam_%s_reajustada" % name,
            lambda forma=name: etrych.family_retuned(forma, ALVOS))
        fora[name + "|intensidade"] = memo(
            "fam_%s_intensidade" % name,
            lambda forma=name: etrych.family_fixed_range(forma, ALVOS))
    fora["lj|c6_sozinho"] = memo("fam_lj_c6_sozinho",
                                 lambda: etrych.sweep_c6_alone(ALVOS))
    return fora


def poco():
    """The solver against the closed forms of the square well.

    Twenty-one depths spread between the first two poles of a, avoiding the
    poles themselves where a diverges. The worst relative deviation of each
    observable is what the notes and the write-up quote.
    """
    pior_a = pior_r0 = 0.0
    for k in range(1, 22):
        v = 0.2 + k * (10.5 - 0.2) / 22.0
        a, r0, nodes = lab.scattering(lab.Well(v, 1.0))
        pior_a = max(pior_a, abs(a / lab.a_well(v, 1.0) - 1.0))
        pior_r0 = max(pior_r0, abs(r0 / lab.r0_well(v, 1.0) - 1.0))
    return {"a": pior_a, "r0": pior_r0}


def ramos():
    """The same (a, r0) on branches carrying 0, 1 and 2 deeply bound dimers.

    Square well only, and the restriction is not cosmetic: the depths are the
    odd squares exactly, so the well is where the branch ladder can be checked
    against a closed form. The qualitative point -- that the branch is a real
    choice with a measurable consequence -- holds for every shape, but the
    quantitative one does NOT: the well keeps nearly the same range on all
    three branches and the other shapes do not, which `alcances_por_ramo`
    measures rather than the caption asserting it.
    """
    row = etrych.TABLE_III[0]
    fora = []
    for deep in (0, 1, 2):
        p1, p2 = etrych.calibrate_on_branch("well", row, deep)
        a, r0, nodes = lab.scattering(lab.Well(p1, p2))
        R = 1.0 / p2
        fora.append({"profundos": deep, "v": p1, "R": R,
                     "R_sobre_rvdw": R / etrych.R_VDW,
                     "a": a, "r0": r0, "estados": nodes})
    return fora


def alcances_por_ramo():
    """How much the range moves between branches, shape by shape.

    The square well barely moves and that is a property of the square well.
    Measured so the write-up can say which it is instead of generalising from
    the one shape its table shows.
    """
    row = etrych.TABLE_III[0]
    fora = {}
    for name in ("well", "mpt", "gauss"):
        Rs = []
        for deep in (0, 1, 2):
            p1, p2 = etrych.calibrate_on_branch(name, row, deep)
            Rs.append(lab.POTENTIALS[name](p1, p2).R)
        fora[name] = {"R": Rs, "espalhamento": max(Rs) / min(Rs) - 1.0}
    return fora


def premissa():
    """What (a, r0) does NOT fix, measured at finite energy and in the dimer.

    The whole construction rests on a slow collision remembering a and r0 and
    nothing else about the potential. That is a statement about finite energy,
    and everything else here is computed at E = 0, so it was assumed and never
    checked. This checks it twice: in the phase shift, and in the one place
    where a bound state exists to be compared.

    The four shapes are matched on (a, r0) at the a_- point and then asked for
    k cot(delta). Each carries a constant offset of its own -- the inner
    boundary is at a finite radius and the cut is not extrapolated here -- so
    the comparison is differenced against a small reference k, and the
    effective-range term is subtracted. What is left is the collision the two
    numbers do not determine.
    """
    # Everything here is returned as LISTS OF PAIRS, never as a dict keyed by
    # a number. This goes through the json cache, and json has only string
    # keys: a fresh run would hand back float keys and a cached run strings,
    # so the formatters below would work once and crash on the second run.
    pots = [lab.POTENTIALS[n](*etrych.at_ratio(n, -0.16)) for n in FORMAS]
    leque = [[k, lab.shape_dependence_at_energy(pots, k)]
             for k in (0.01, 0.02, 0.05, 0.10)]

    # the leftover of one shape, and the power it follows: if it were an
    # artefact of the matching it would not shrink with k at all
    pot = lab.Well(*etrych.at_ratio("well", -0.16))
    a, r0, nodes = lab.scattering(pot)
    base = lab.k_cot_delta(pot, 1e-3)
    sobras = []
    for k in (0.01, 0.02, 0.05):
        medido = lab.k_cot_delta(pot, k) - base
        sobras.append([k, abs(medido - 0.5 * r0 * (k * k - 1e-6))])
    potencia = math.log(sobras[-1][1] / sobras[0][1]) / math.log(5.0)
    coef = sobras[-1][1] / sobras[-1][0] ** 4

    # and the bound sector: the shallow dimer, where a > 0 puts one
    dimeros = []
    for alvo in (0.06, 0.16):
        Es = []
        for name in FORMAS:
            p = lab.POTENTIALS[name](*etrych.at_ratio(name, alvo))
            n, energias = lab.count_bound_states(p)
            Es.append(energias[0])
        # the energies are negative, so max/min would come back below one and
        # the spread negative. Compare binding energies, which are magnitudes.
        fundo = [abs(E) for E in Es]
        dimeros.append([alvo, {"E": Es,
                               "espalhamento": max(fundo) / min(fundo) - 1.0}])

    return {"leque": leque, "potencia": potencia, "coeficiente": coef,
            "sobras": sobras, "dimeros": dimeros}


def escadas():
    """The branch boundaries against the two closed forms they have.

    The boundary is the pole of a, located by bisection. For two of the four
    shapes the ladder of poles is analytic, which turns the search into
    something checkable rather than merely convergent.
    """
    fora = {}
    for name, formula, rotulo in (
            ("well", lambda n: float((2 * n - 1) ** 2),
             r"\ensuremath{(2n-1)^2}"),
            ("mpt", lambda n: float(n * (2 * n - 1)),
             r"\ensuremath{n(2n-1)}")):
        f = etrych.fronteiras(name, teto=60.0)
        passos = etrych.degraus(name, teto=60.0)
        fora[name] = {
            "rotulo": rotulo,
            "polos": [f[n] for n in sorted(f)],
            "exatos": [formula(n) for n in sorted(f)],
            "erro": max(abs(f[n] / formula(n) - 1.0) for n in f),
            "atraso": passos[min(passos)] - f[min(f)],
        }
    return fora


def positivas():
    """The calibration at a = a_p, and the shape dependence there.

    Table III publishes the positive-side feature for four of its five rows,
    and the experiment reports it as the easier of the two to measure. Whether
    that helps depends on the theory side being no worse conditioned there,
    which is what the spread between shapes says.
    """
    linhas = []
    for row in etrych.TABLE_III:
        r = etrych.ratio_p(row)
        if r is None:
            continue
        for name in FORMAS:
            p1, p2 = etrych.calibrate_p(name, row)
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            linhas.append({"B_res": row["B_txt"], "r_eff": row["r_eff"],
                           "a_p": row["a_p"], "razao": r, "forma": name,
                           "p1": p1, "p2": p2, "a": a, "r0": r0,
                           "estados": nodes})
    leque = {}
    for row in etrych.TABLE_III:
        r = etrych.ratio_p(row)
        if r is None:
            continue
        v = [etrych.family_fixed_range(n, [r])[0]["r0"] for n in FORMAS]
        leque[row["B_txt"]] = max(v) / min(v) - 1.0
    return {"linhas": linhas, "leque": leque}


def tabela_ramos(rm):
    return ["%d & %.5f & %.1f & %.2f & $%.3f$ & %.3f \\\\"
            % (l["profundos"], l["v"], l["R"], l["R_sobre_rvdw"], l["a"], l["r0"])
            for l in rm]


def tabela_positivas(pos):
    """One row per shape per published a_p.

    The parameters go through notacao() for the same reason the a_- table does:
    the Lennard-Jones coefficients are 1e8 and 1e18 at these ranges, and a
    fixed-point format turns them into unreadable walls of digits.
    """
    linhas = []
    visto = None
    for l in pos["linhas"]:
        if l["B_res"] != visto:
            if visto is not None:
                linhas.append(r"\addlinespace")
            cabeca = "%s & $%d$" % (l["B_res"], l["a_p"])
        else:
            cabeca = " & "
        visto = l["B_res"]
        linhas.append("%s & %s & %s & %s & $%.3f$ & %.3f & %d \\\\"
                      % (cabeca, MACRO[l["forma"]], notacao(l["p1"]),
                         notacao(l["p2"]), l["a"], l["r0"], l["estados"]))
    return linhas


def tabela_premissa(pr):
    """The part of the collision (a, r0) does not fix, against k r0."""
    return ["$%.2f$ & %s \\\\" % (k, notacao(v)) for k, v in pr["leque"]]


def cauda():
    """r0 of one fixed Lennard-Jones as the tail cut is pushed outwards.

    The other three shapes have exponential or faster tails and stop moving at
    once. This is where the error of this calculation lives, so it is reported
    rather than hidden.
    """
    p1, p2 = etrych.unit_solution("lj")
    pot = lab.LennardJones(p1, p2)
    base = pot.R
    fora = []
    for m in (1, 3, 9, 27):
        a, r0, nodes = lab._num_corte(pot, base * m)
        fora.append({"corte": "%dR" % m, "R": base * m,
                     "inv_a": 1.0 / a, "r0": r0})
    return fora


def residuo_da_cauda():
    """How far the tail extrapolation is from itself, one factor out.

    Extrapolating from (R, 3R) and from (3R, 9R) shares no cut, so the
    difference between them is the honest error the tail treatment leaves --
    not an estimate of it. Measured at the a_- point, where the calibration
    actually sits, rather than at unitarity.

    Also returns the (9R, 27R) extrapolation, which is WORSE: past a few
    hundred times the range the asymptotic line, growing as r, carries the
    accumulated round-off and r0 starts to collapse. So 3R is not laziness, it
    is where the two error sources cross.
    """
    pot = lab.LennardJones(*etrych.at_ratio("lj", -0.16))
    def extrap(f):
        perto = lab._num_corte(pot, f * pot.R)[1]
        longe = lab._num_corte(pot, 3.0 * f * pot.R)[1]
        return (3.0 * longe - perto) / 2.0
    um, tres, nove = extrap(1.0), extrap(3.0), extrap(9.0)
    return {"r0": um,
            "residuo": abs(tres / um - 1.0),
            "residuo_longe": abs(nove / um - 1.0)}


# =============================================================================
#  What gets written
# =============================================================================
def tabela_iii():
    """The transcription itself, so the notes and the code cannot disagree."""
    linhas = []
    for row in etrych.TABLE_III:
        mais = (r"%d & $\pm$ & %d" % (row["a_p"], row["a_p_err"])
                if row["a_p"] else r"\multicolumn{3}{c}{---}")
        linhas.append(
            r"$|%s\rangle$ & %s & %.1f & %d & $%d$ & $\pm$ & %d & %s & $%.4f$ \\"
            % (row["state"].strip("|>").replace(",", ","), row["B_txt"],
               row["s_res"], row["r_eff"], row["a_minus"], row["a_minus_err"],
               mais, etrych.ratio(row)))
    return linhas


def tabela_unitaria(tabela):
    linhas = []
    for linha in tabela:
        if linha["calibracao"] != "unitaria" or linha["r_eff"] != 136:
            continue
        p1, p2 = etrych.unit_solution(linha["forma"])
        a, r0, n = lab.scattering(lab.POTENTIALS[linha["forma"]](p1, p2))
        linhas.append("%s & %.9f & %.9f & %s & %.9f \\\\"
                      % (MACRO[linha["forma"]], p1, p2, notacao(1.0 / a), r0))
    return linhas


def tabela_aminus(tabela):
    linhas = []
    for row in etrych.TABLE_III:
        primeiro = True
        for linha in tabela:
            if linha["calibracao"] != "a_menos" or linha["r_eff"] != row["r_eff"] \
                    or linha["a_menos"] != row["a_minus"]:
                continue
            cabeca = ("%d & $%d$" % (row["r_eff"], row["a_minus"])) if primeiro else " & "
            primeiro = False
            linhas.append("%s & %s & %s & %s & $%.3f$ & %.3f & %d \\\\"
                          % (cabeca, MACRO[linha["forma"]], notacao(linha["p1"]),
                             notacao(linha["p2"]), linha["a"], linha["r0"],
                             linha["estados"]))
        linhas.append(r"\addlinespace")
    return linhas


def tabela_familias(fam):
    """Per shape: how far r0 drifts along the walk, and how far the parameter
    has to move. NOT the shape dependence -- that is the spread BETWEEN shapes
    at the same point, which is a different number and is emitted separately
    below. The two were confused once; both are computed here so neither has
    to be recalled.
    """
    linhas, deriva, forca = [], {}, {}
    for name in FORMAS:
        fixa = fam[name + "|intensidade"]
        r0f = [l["r0"] for l in fixa]
        raz = [l["strength_ratio"] for l in fixa]
        r0r = [l["r0"] for l in fam[name + "|reajustada"]]
        deriva[name] = max(r0f) / min(r0f) - 1.0
        forca[name] = max(raz) / min(raz) - 1.0
        linhas.append("%s & %.6f & %.6f & %.1f\\%% & %.8f & %.8f & %.1f\\%% \\\\"
                      % (MACRO[name], min(r0f), max(r0f), 100.0 * deriva[name],
                         min(r0r), max(r0r), 100.0 * forca[name]))
    c6 = fam["lj|c6_sozinho"]
    r0e = [l["r0"] for l in c6]
    raze = [l["strength_ratio"] for l in c6]
    rme = [l["r_min"] for l in c6]
    deriva["c6"] = max(r0e) / min(r0e) - 1.0
    forca["c6"] = max(raze) / min(raze) - 1.0
    deriva["rmin"] = max(rme) / min(rme) - 1.0
    linhas.append(r"\addlinespace")
    linhas.append(r"\quad {\itshape \FormaCseis} & %.6f & %.6f & %.1f\%% &"
                  r" \multicolumn{2}{c}{---} & %.1f\%% \\"
                  % (min(r0e), max(r0e), 100.0 * deriva["c6"], 100.0 * forca["c6"]))
    return linhas, deriva, forca


def leque(fam):
    """The spread BETWEEN the four shapes at the same point: the shape
    dependence proper. Zero at unitarity by construction -- all four are
    calibrated there -- and growing close to linearly away from it.
    """
    fora = {}
    for k, inv in enumerate(ALVOS):
        v = [fam[name + "|intensidade"][k]["r0"] for name in FORMAS]
        fora[round(inv, 10)] = max(v) / min(v) - 1.0
    return fora


def macros(tabela, deriva, forca, lq):
    v, mu = etrych.unit_solution("well")
    exato = math.pi ** 2 / 8.0
    unit = [l for l in tabela if l["calibracao"] == "unitaria"]
    alvo = [l for l in tabela if l["calibracao"] == "a_menos"]
    linhas = [
        r"\newcommand{\EtrychVunit}{%.10f}" % v,
        r"\newcommand{\EtrychVexato}{%.10f}" % exato,
        r"\newcommand{\EtrychVerro}{%s}" % notacao(abs(v / exato - 1.0)),
        r"\newcommand{\EtrychUnitInv}{%s}"
        % notacao(max(abs(l["inv_a"]) for l in unit)),
        r"\newcommand{\EtrychUnitRzero}{%s}"
        % notacao(max(abs(l["r0"] / l["r_eff"] - 1.0) for l in unit)),
        r"\newcommand{\EtrychAlvoA}{%s}"
        % notacao(max(abs(l["a"] / l["a_menos"] - 1.0) for l in alvo)),
        r"\newcommand{\EtrychAlvoRzero}{%s}"
        % notacao(max(abs(l["r0"] / l["r_eff"] - 1.0) for l in alvo)),
        # the ratio the premise is measured at, so the caption does not type it
        r"\newcommand{\EtrychAlvoRazao}{%.2f}" % -0.16,
        r"\newcommand{\EtrychAlvoBaixo}{%.2f}" % ALVOS[0],
        r"\newcommand{\EtrychAlvoAlto}{%.2f}" % ALVOS[-1],
        r"\newcommand{\EtrychLequeBorda}{%.1f\%%}" % (100.0 * lq[ALVOS[0]]),
        r"\newcommand{\EtrychLequeAmenos}{%.1f\%%}" % (100.0 * lq[-0.16]),
        r"\newcommand{\EtrychLequeRaso}{%.1f\%%}" % (100.0 * lq[-0.06]),
        r"\newcommand{\EtrychMinimoCseis}{%.1f\%%}" % (100.0 * deriva["rmin"]),
    ]
    for name in FORMAS:
        chave = name.capitalize() if name != "mpt" else "Mpt"
        linhas.append(r"\newcommand{\EtrychDeriva%s}{%.1f\%%}"
                      % (chave, 100.0 * deriva[name]))
        linhas.append(r"\newcommand{\EtrychForca%s}{%.1f\%%}"
                      % (chave, 100.0 * forca[name]))
    linhas.append(r"\newcommand{\EtrychDerivaCseis}{%.1f\%%}" % (100.0 * deriva["c6"]))
    linhas.append(r"\newcommand{\EtrychForcaCseis}{%.1f\%%}" % (100.0 * forca["c6"]))
    return linhas


def escreve_tex():
    caminho = os.path.join(AQUI, "tabelas", "gerado.tex")
    if not os.path.isdir(os.path.dirname(caminho)):
        os.makedirs(os.path.dirname(caminho))
    with open(caminho, "w", encoding="utf-8") as f:
        f.write("%% generated by gerar.py -- do not edit\n")
        for macro, linhas in SAIDA_TEX:
            if macro is None:
                f.write("\n".join(linhas) + "\n")
            else:
                f.write("\n\\newcommand{\\%s}{%%\n" % macro)
                f.write("\n".join(linhas) + "\n}\n")
    print("  wrote", caminho)


def escreve_csv(nome, linhas):
    caminho = os.path.join(SAIDA, nome)
    # lineterminator: the csv module writes \r\n by default, which lands in a
    # repository that is otherwise all LF.
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()),
                           lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    print("  wrote", caminho)


def figuras(fam):
    pasta = os.path.join(SAIDA, "figuras")
    if not os.path.isdir(pasta):
        os.makedirs(pasta)
    razoes = [etrych.ratio(row) for row in etrych.TABLE_III]
    positivas = [r for r in (etrych.ratio_p(row) for row in etrych.TABLE_III)
                 if r is not None]

    def marcas(ax):
        for x in razoes:
            ax.axvline(x, color="#566573", ls=":", lw=0.7)
        for x in positivas:
            ax.axvline(x, color="#b7950b", ls=":", lw=0.7)
        ax.axvline(0.0, color="k", lw=0.5)
        ax.set_xlabel(r"$r_{\rm eff}/a$")
        ax.grid(alpha=0.2)

    fig, ax = plt.subplots(figsize=(6.2, 4.0), layout="tight")
    for name in FORMAS:
        ax.plot(ALVOS, [l["r0"] for l in fam[name + "|reajustada"]],
                "--", color=COR[name], lw=1.2)
        ax.plot(ALVOS, [l["r0"] for l in fam[name + "|intensidade"]],
                "-", color=COR[name], lw=1.6, label=MACRO[name][6:])
    ax.plot(ALVOS, [l["r0"] for l in fam["lj|c6_sozinho"]],
            ":", color=COR["lj"], lw=1.6, label="Lj, $C_6$ so")
    marcas(ax)
    ax.set(ylabel=r"$r_0/r_{\rm eff}$",
           title="tracejado: ambos reajustados   cheio: so a intensidade")
    ax.legend(fontsize=7.5, loc="upper right")
    fig.savefig(os.path.join(pasta, "alcance_efetivo.pdf"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0), layout="tight")
    for name in FORMAS:
        ax.plot(ALVOS, [l["strength_ratio"] for l in fam[name + "|intensidade"]],
                "-", color=COR[name], lw=1.6, label=MACRO[name][6:])
    ax.plot(ALVOS, [l["strength_ratio"] for l in fam["lj|c6_sozinho"]],
            ":", color=COR["lj"], lw=1.6, label="Lj, $C_6$ so")
    marcas(ax)
    ax.set(ylabel="intensidade / intensidade na unitariedade",
           title="o parametro ao longo do caminho")
    ax.legend(fontsize=7.5, loc="upper left")
    fig.savefig(os.path.join(pasta, "intensidade.pdf"))
    plt.close(fig)

    # The four shapes give the same count at every point, so drawn plainly only
    # the last is visible -- which reads as though the others were not
    # computed. Staggered widths make the agreement something the figure shows.
    fig, ax = plt.subplots(figsize=(6.2, 3.2), layout="tight")
    for k, name in enumerate(FORMAS):
        ax.step(ALVOS, [l["nodes"] for l in fam[name + "|intensidade"]],
                where="post", color=COR[name], lw=5.0 - 1.1 * k, alpha=0.85,
                label=MACRO[name][6:])
    marcas(ax)
    ax.set(ylabel="estados ligados", yticks=[0, 1, 2],
           title="um estado entra na ressonancia")
    ax.legend(fontsize=7.5, loc="upper left")
    fig.savefig(os.path.join(pasta, "estados_ligados.pdf"))
    plt.close(fig)
    print("  wrote", pasta)


def main():
    if not os.path.isdir(SAIDA):
        os.makedirs(SAIDA)

    print("tables")
    tabela = memo("calibracoes", calibracoes)
    fam = familias()
    linhas_fam, deriva, forca = tabela_familias(fam)
    lq = leque(fam)

    SAIDA_TEX.append(("TabEtrychIIICorpo", tabela_iii()))
    SAIDA_TEX.append(("TabEtrychUnitarioCorpo", tabela_unitaria(tabela)))
    SAIDA_TEX.append(("TabEtrychAminusCorpo", tabela_aminus(tabela)))
    SAIDA_TEX.append(("TabEtrychFamiliasCorpo", linhas_fam))
    pc = memo("poco", poco)
    rm = memo("ramos", ramos)
    SAIDA_TEX.append(("TabEtrychRamosCorpo", tabela_ramos(rm)))

    pr = memo("premissa", premissa)
    esc = memo("escadas", escadas)
    pos = memo("positivas", positivas)
    alc = memo("alcances_ramo", alcances_por_ramo)
    rc = memo("residuo_cauda", residuo_da_cauda)
    SAIDA_TEX.append(("TabEtrychPositivasCorpo", tabela_positivas(pos)))
    SAIDA_TEX.append(("TabEtrychPremissaCorpo", tabela_premissa(pr)))

    etas = [row["eta_minus"] for row in etrych.TABLE_III]
    extra = [r"\newcommand{\EtrychEtaMenor}{%.2f}" % min(etas),
             r"\newcommand{\EtrychEtaMaior}{%.2f}" % max(etas),
             r"\newcommand{\EtrychRamoAlcanceBaixo}{%.2f}"
             % min(l["R_sobre_rvdw"] for l in rm),
             r"\newcommand{\EtrychRamoAlcanceAlto}{%.2f}"
             % max(l["R_sobre_rvdw"] for l in rm),
             r"\newcommand{\EtrychPocoA}{%s}" % notacao(pc["a"]),
             r"\newcommand{\EtrychPocoRzero}{%s}" % notacao(pc["r0"]),
             # the premise
             r"\newcommand{\EtrychPremissaBaixo}{%s}" % notacao(pr["leque"][0][1]),
             r"\newcommand{\EtrychPremissaMeio}{%s}" % notacao(pr["leque"][2][1]),
             r"\newcommand{\EtrychPremissaAlto}{%s}" % notacao(pr["leque"][3][1]),
             r"\newcommand{\EtrychPremissaKbaixo}{%.2f}" % pr["leque"][0][0],
             r"\newcommand{\EtrychPremissaKmeio}{%.2f}" % pr["leque"][2][0],
             r"\newcommand{\EtrychPremissaKalto}{%.2f}" % pr["leque"][3][0],
             r"\newcommand{\EtrychPremissaPotencia}{%.3f}" % pr["potencia"],
             r"\newcommand{\EtrychPremissaCoef}{%.3f}" % pr["coeficiente"],
             r"\newcommand{\EtrychDimeroLeque}{%s}"
             % notacao(pr["dimeros"][0][1]["espalhamento"]),
             r"\newcommand{\EtrychDimeroRazao}{%.2f}" % pr["dimeros"][0][0],
             r"\newcommand{\EtrychDimeroEnergia}{%s}"
             % notacao(abs(pr["dimeros"][0][1]["E"][0])),
             # the ladders
             r"\newcommand{\EtrychEscadaPoco}{%s}" % esc["well"]["rotulo"],
             r"\newcommand{\EtrychEscadaMpt}{%s}" % esc["mpt"]["rotulo"],
             r"\newcommand{\EtrychEscadaErro}{%s}"
             % notacao(max(esc[n]["erro"] for n in esc)),
             r"\newcommand{\EtrychEscadaAtraso}{%s}" % notacao(esc["well"]["atraso"]),
             # the positive side
             r"\newcommand{\EtrychPositivasN}{%d}" % len(pos["linhas"]),
             r"\newcommand{\EtrychPositivaA}{%s}"
             % notacao(max(abs(l["a"] / l["a_p"] - 1.0) for l in pos["linhas"])),
             r"\newcommand{\EtrychPositivaRzero}{%s}"
             % notacao(max(abs(l["r0"] / l["r_eff"] - 1.0) for l in pos["linhas"])),
             r"\newcommand{\EtrychLequePositivoBaixo}{%.1f\%%}"
             % (100.0 * min(pos["leque"].values())),
             r"\newcommand{\EtrychLequePositivoAlto}{%.1f\%%}"
             % (100.0 * max(pos["leque"].values())),
             r"\newcommand{\EtrychLequeMaisDezesseis}{%.1f\%%}" % (100.0 * lq[0.16]),
             # how much the range moves between branches, by shape
             r"\newcommand{\EtrychRamoPocoMove}{%.0f\%%}"
             % (100.0 * alc["well"]["espalhamento"]),
             r"\newcommand{\EtrychRamoOutrasMove}{%.1f}"
             % max(1.0 + alc[n]["espalhamento"] for n in ("mpt", "gauss")),
             # the honest error the tail treatment leaves
             r"\newcommand{\EtrychCaudaResiduo}{%s}" % notacao(rc["residuo"]),
             r"\newcommand{\EtrychCaudaResiduoLonge}{%s}"
             % notacao(rc["residuo_longe"])]
    SAIDA_TEX.append((None, macros(tabela, deriva, forca, lq) + extra))
    escreve_tex()

    print("\ndata")
    escreve_csv("calibracoes.csv", tabela)
    plano = []
    for chave, registros in fam.items():
        name, prescricao = chave.split("|")
        for reg in registros:
            plano.append({"prescricao": prescricao, "forma": name,
                          "r_eff_sobre_a": reg["inv_a"],
                          "p1": reg["p1"], "p2": reg["p2"],
                          "a": reg["a"], "r0": reg["r0"],
                          "estados": reg["nodes"],
                          "fator_intensidade": reg.get("strength_ratio", ""),
                          "r_minimo": reg.get("r_min", "")})
    escreve_csv("familias.csv", plano)
    escreve_csv("corte_da_cauda.csv", memo("cauda", cauda))

    print("\nfigures")
    figuras(fam)

    print("\nthe shape dependence, between shapes at the same point")
    for inv in (ALVOS[0], -0.16, -0.06, 0.0, 0.16, ALVOS[-1]):
        print("  r_eff/a = %+5.2f   %5.2f%%" % (inv, 100.0 * lq[round(inv, 10)]))


if __name__ == "__main__":
    main()

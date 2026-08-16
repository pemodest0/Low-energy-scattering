# -*- coding: utf-8 -*-
"""
Gera os fragmentos LaTeX (tabelas) do relatório a partir dos CSVs em
resultados/, copia as figuras para relatorio/figuras/ e compila o PDF.
"""
import csv
import os
import shutil
import subprocess

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_RES = os.path.join(AQUI, "resultados")
DIR_FIG = os.path.join(AQUI, "figuras")
DIR_REL = os.path.join(AQUI, "relatorio")

CASOS = {"nn": "n--n", "unitario": "unit.", "deuteron": "d\\^euteron"}
POTS = {"poco": "po\\c{c}o", "mpt": "mPT", "gauss": "gaussiano",
        "lj": "Lennard-Jones"}


def _le(nome):
    with open(os.path.join(DIR_RES, nome), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _fmt(x, casas=4):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "---"
    if v != 0 and (abs(v) >= 1e4 or abs(v) < 1e-3):
        mant, exp = f"{v:.2e}".split("e")
        return f"${mant} \\times 10^{{{int(exp)}}}$"
    return f"${v:.{casas}f}$"


def tabela_validacao():
    linhas = [r"\begin{tabular}{llcccc}", r"\hline",
              r"potencial & caso & método & $a$ num.\ (fm) & $a$ anal.\ (fm) & erro rel. \\",
              r"\hline"]
    for x in _le("validacao_analitica.csv"):
        linhas.append(
            f"{POTS[x['potencial']]} & {CASOS[x['caso']]} & {x['metodo']} & "
            f"{_fmt(x['a_num'])} & {_fmt(x['a_analitico'])} & "
            f"{_fmt(x['erro_rel_a'])} \\\\")
    linhas += [r"\hline", r"\end{tabular}"]
    return "\n".join(linhas)


def tabela_reproducao():
    linhas = [r"\begin{tabular}{llccccc}", r"\hline",
              r"caso & potencial & $a$ num.\ (fm) & $a$ artigo & "
              r"$r_0$ num.\ (fm) & $r_0$ artigo & nós \\", r"\hline"]
    for x in _le("reproducao_tabelas_3_4.csv"):
        linhas.append(
            f"{CASOS[x['caso']]} & {POTS[x['potencial']]} & "
            f"{_fmt(x['a_num'])} & {_fmt(x['a_artigo'])} & "
            f"{_fmt(x['r0_num'])} & {_fmt(x['r0_artigo'], 2)} & "
            f"{x['nos']} ({x['nos_esperado']}) \\\\")
    linhas += [r"\hline", r"\end{tabular}"]
    return "\n".join(linhas)


def tabela_ajustes():
    linhas = [r"\begin{tabular}{llcccccc}", r"\hline",
              r"caso & potencial & $p_1$ ajust. & $p_1$ artigo & "
              r"$p_2$ ajust. & $p_2$ artigo & $a$ (fm) & $r_0$ (fm) \\",
              r"\hline"]
    for x in _le("ajustes.csv"):
        pot = x["potencial"]
        n1 = "C6" if pot == "lj" else "v"
        n2 = "C12" if pot == "lj" else "mu"
        linhas.append(
            f"{CASOS[x['caso']]} & {POTS[pot]} & "
            f"{_fmt(x[f'{n1}_ajustado'], 5)} & {_fmt(x[f'{n1}_artigo'], 5)} & "
            f"{_fmt(x[f'{n2}_ajustado'], 5)} & {_fmt(x[f'{n2}_artigo'], 5)} & "
            f"{_fmt(x['a_obtido'])} & {_fmt(x['r0_obtido'])} \\\\")
    linhas += [r"\hline", r"\end{tabular}"]
    return "\n".join(linhas)


def tabela_energias():
    linhas = [r"\begin{tabular}{lccccc}", r"\hline",
              r"sistema & $E_{\rm zr}$ & $E_{\rm zr}$ artigo & "
              r"$E_{\rm fr}$ & $E_{\rm fr}$ artigo & $E$ exp. \\", r"\hline"]
    rot = {"dimero_4He": r"d\'imero $^4$He (mK)", "deuteron": r"d\^euteron (MeV)"}
    for x in _le("energias_tabela1.csv"):
        fator = 1e3 if x["sistema"] == "dimero_4He" else 1.0
        linhas.append(
            f"{rot[x['sistema']]} & {_fmt(float(x['E_zr']) * fator, 3)} & "
            f"{_fmt(float(x['E_zr_artigo']) * fator, 3)} & "
            f"{_fmt(float(x['E_fr']) * fator, 3)} & "
            f"{_fmt(float(x['E_fr_artigo']) * fator, 3)} & "
            f"{_fmt(float(x['E_ref']) * fator, 3)} \\\\")
    linhas += [r"\hline", r"\end{tabular}"]
    return "\n".join(linhas)


def gerar_fragmentos():
    frag = os.path.join(DIR_REL, "tabelas_geradas.tex")
    with open(frag, "w", encoding="utf-8") as f:
        f.write("% Arquivo gerado automaticamente por src/gerar_relatorio.py\n")
        f.write("\\newcommand{\\tabvalidacao}{%\n" + tabela_validacao() + "}\n")
        f.write("\\newcommand{\\tabreproducao}{%\n" + tabela_reproducao() + "}\n")
        f.write("\\newcommand{\\tabajustes}{%\n" + tabela_ajustes() + "}\n")
        f.write("\\newcommand{\\tabenergias}{%\n" + tabela_energias() + "}\n")
    print("  ->", frag)


def copiar_figuras():
    destino = os.path.join(DIR_REL, "figuras")
    os.makedirs(destino, exist_ok=True)
    for nome in os.listdir(DIR_FIG):
        if nome.endswith(".png"):
            shutil.copy(os.path.join(DIR_FIG, nome), destino)
    print("  -> figuras copiadas para", destino)


def compilar():
    for _ in range(2):                       # duas passadas (referências)
        subprocess.run(["pdflatex", "-interaction=nonstopmode",
                        "relatorio.tex"],
                       cwd=DIR_REL, capture_output=True, check=False)
    pdf = os.path.join(DIR_REL, "relatorio.pdf")
    if os.path.exists(pdf):
        print("  ->", pdf)
    else:
        raise RuntimeError("pdflatex não gerou o PDF; veja relatorio.log")


def gerar(compila=True):
    gerar_fragmentos()
    copiar_figuras()
    if compila:
        compilar()


if __name__ == "__main__":
    gerar()

# -*- coding: utf-8 -*-
"""
Pipeline ponta-a-ponta do laboratório de espalhamento.

Uso:
    python main.py                # tudo: resultados + figuras + relatório
    python main.py --sem-ajustes  # pula o ajuste de parâmetros (etapa lenta)
    python main.py --so-figuras   # apenas as figuras
    python main.py --so-relatorio # apenas recompila o relatório

Saídas:
    resultados/  CSVs + resumo.md
    figuras/     figs 7-10 + convergência (.png)
    relatorio/   relatorio.pdf

Tempo estimado (laptop comum): ~1 min sem ajustes; os ajustes completos
(incluindo os 3 casos do Lennard-Jones) adicionam ~3-5 min.
"""
import argparse
import sys
import time


def main():
    ap = argparse.ArgumentParser(description="Laboratório de espalhamento (RBEF 45, e20230079)")
    ap.add_argument("--sem-ajustes", action="store_true",
                    help="pula a etapa (lenta) de ajuste de parâmetros")
    ap.add_argument("--so-figuras", action="store_true")
    ap.add_argument("--so-relatorio", action="store_true")
    args = ap.parse_args()

    t0 = time.time()

    if not (args.so_figuras or args.so_relatorio):
        from src import gerar_resultados as gr
        print("== etapa 1: resultados numéricos ==", flush=True)
        va = gr.validacao_analitica(1e-3)
        gr._escreve_csv("validacao_analitica.csv", va)
        rt = gr.reproducao_tabelas(1e-3)
        gr._escreve_csv("reproducao_tabelas_3_4.csv", rt)
        e1 = gr.energias_tabela1()
        gr._escreve_csv("energias_tabela1.csv", e1)
        cv = gr.convergencia()
        gr._escreve_csv("convergencia.csv", cv)
        if args.sem_ajustes:
            print("(ajustes pulados; usando ajustes.csv existente, se houver)")
            aj = []
        else:
            print("== etapa 2: ajustes de parâmetros (lento) ==", flush=True)
            aj = gr.rodar_ajustes(dr=1e-3)
            gr._escreve_csv("ajustes.csv", aj)
        gr._escreve_resumo(va, rt, e1, aj)

    if not args.so_relatorio:
        from src import gerar_figuras as gf
        print("== etapa 3: figuras ==", flush=True)
        gf.gerar_todas()

    if not args.so_figuras:
        from src import gerar_relatorio as grl
        print("== etapa 4: relatório LaTeX -> PDF ==", flush=True)
        try:
            grl.gerar(compila=True)
        except (RuntimeError, FileNotFoundError) as exc:
            print(f"[aviso] compilação falhou ({exc}); o .tex está pronto em relatorio/")

    print(f"concluído em {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())

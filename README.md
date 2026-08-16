# Laboratório numérico de espalhamento de baixa energia

Rascunho **v0** — reprodução integral de:

> M. Macêdo-Lima e L. Madeira, *Scattering length and effective range of
> microscopic two-body potentials*, Rev. Bras. Ens. Fís. **45**, e20230079
> (2023). DOI: [10.1590/1806-9126-RBEF-2023-0079](https://doi.org/10.1590/1806-9126-RBEF-2023-0079)

Calcula o comprimento de espalhamento `a` e o alcance efetivo `r0` de
potenciais de dois corpos (onda-s, energia zero), com dois métodos
numéricos, validação analítica, ajuste de parâmetros, figuras, laboratório
interativo e relatório em PDF.

## Dependências

Python ≥ 3.9 com `numpy`, `scipy` e `matplotlib` (v1; a v0 era sem scipy). Para o relatório:
`pdflatex` no PATH (opcional — sem ele o `.tex` fica pronto).

## Roadmap do projeto

Abra **`ROADMAP.html`** no navegador: linha do tempo interativa com de onde
partimos, onde estamos, onde queremos chegar, o que deu certo, o que deu
errado (e a lição), o que precisa melhorar, e um checklist de próximos
passos que fica salvo no navegador.

## Por onde começar (se está chegando agora)

**Dois cliques em `app.bat`** — o ambiente unificado no navegador: 7
estações (espalhamento → ajuste → Schrödinger → armadilha 3D→2D → Feshbach
³⁹K → Efimov → validação), cada uma com as abas Simular / Teoria / Código /
Resultados e tempo estimado antes de rodar.

Para leitura: `INICIO.html` (hub com tudo embutido). Para estudo:
`estudar.bat` — **comece pelo notebook `00_do_zero`** (não assume nada,
vai devagar), depois 01–04. Sliders soltos: `laboratorio.bat`. Todos
instalam o que faltar sozinhos na primeira vez.

## Diário de bordo

`diario/DIARIO.md` registra, por dia: mudanças, descobertas, erros
encontrados (com a lição) e a versão. Nova entrada: `python -m src.diario`.

## Instalação

```bash
pip install -r requirements.txt      # ou: pip install -e .
```

## Como rodar tudo

```bash
python main.py                 # resultados + ajustes + figuras + PDF (~5 min)
python main.py --sem-ajustes   # versão rápida (~1 min)
python interativo/laboratorio.py   # laboratório com sliders
```

## Estrutura

```
laboratorio_espalhamento/
├── main.py                  # pipeline ponta-a-ponta
├── src/
│   ├── constantes.py        # hbar c, massas (validação da Tabela 1)
│   ├── potenciais.py        # poço, mPT, gaussiano, Lennard-Jones
│   ├── solvers.py           # diferença central (Eq. 99) e Numerov (Eq. 101)
│   ├── espalhamento.py      # a (Eq. 110), r0 (Eq. 56, trapézio+Simpson), nós
│   ├── analitico.py         # Eq. 80/92 (poço), Eq. 117/118 (mPT), E_zr/E_fr
│   ├── ajuste.py            # laços aninhados (Seção 4.5), Illinois
│   ├── tabelas_artigo.py    # alvos e referências das Tabelas 1-4
│   ├── gerar_resultados.py  # CSVs + resumo.md
│   ├── gerar_figuras.py     # Figs. 7-10 + convergência
│   └── gerar_relatorio.py   # tabelas LaTeX + compilação do PDF
├── figuras/                 # .png gerados
├── resultados/              # .csv + resumo.md
├── notebooks/               # visita guiada executável
├── interativo/laboratorio.py
├── relatorio/               # relatorio.tex + relatorio.pdf
└── README.md
```

## Convenções físicas (importante)

- Unidades adimensionais: `hbar = m_r = 1`, comprimentos em fm
  (`epsilon = hbar²/(m_r l²)`, `l = 1 fm`). Equação resolvida:
  `u'' = 2 V̄ u` com `u(0)=0`, `u(Δr)=1`.
- Alcance numérico: `|V̄(R)| ≤ 1e-15`; LJ começa em `r_min` com
  `U(r_min) ≈ 1e10` e `u=0` no caroço.
- **Nota (fator 2 no LJ):** as constantes da Tabela 4 do artigo só
  reproduzem os `(a, r0)` publicados com `u'' = (C12/r¹² − C6/r⁶) u`,
  i.e. `V_LJ = (hbar²/2m_r)[...]` — e não `hbar²/m_r` como impresso na
  Eq. (121). Ver docstring em `src/potenciais.py` e o relatório.
- **Nota (Numerov × descontinuidade):** para o poço esférico a
  descontinuidade na borda degrada o Numerov para erro O(Δr); a diferença
  central com a borda tratada como ponto interior fica O(Δr²). Para
  potenciais suaves o Numerov é ordens de grandeza melhor (ver
  `figuras/fig_convergencia.png`).

## Critérios de aceite (status)

- [x] Dois métodos concordam entre si e com o analítico (poço: ~3e-7;
      mPT: ~1e-10 com Numerov; ver `resultados/validacao_analitica.csv`)
- [x] Parâmetros das Tabelas 3 e 4 reproduzidos (~0,1%; LJ ~1% com
      degenerescência do par (C12, C6); ver `resultados/ajustes.csv`)
- [x] Figuras 7–10 equivalentes às do artigo (`figuras/`)
- [x] Energias da Tabela 1: E_zr = −1,48 mK / −1,416 MeV;
      E_fr = −1,63 mK / −2,223 MeV (`resultados/energias_tabela1.csv`)
- [x] `main.py` roda tudo ponta-a-ponta

Este é um rascunho v0: a versão final será reescrita linha a linha.

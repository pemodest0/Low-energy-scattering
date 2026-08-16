# MAPA — onde está cada coisa

Uma página só. Se você abrir a pasta e não souber por onde começar, é aqui.

**Regra da organização:** tudo é agrupado por **número de corpos** — um corpo,
dois corpos, três corpos — porque é assim que a física vai ficando mais difícil,
e é assim que a dissertação vai andar.

---

## O caminho de aprendizado (na ordem)

| Passo | Abra isto | O que você aprende |
|---|---|---|
| 1 | `notebooks/1_um_corpo/04_schrodinger_canonica.ipynb` | Uma partícula num potencial: oscilador, hidrogênio |
| 2 | `notebooks/2_dois_corpos/LABORATORIO_2CORPOS_COMPLETO.ipynb` | **O principal.** 21 exercícios, você escreve cada linha |
| 3 | `notas_teoria/DoZeroAoTrimero_ParteI.pdf` | A teoria de dois corpos por escrito |
| 4 | `notas_teoria/DoZeroAoTrimero_ParteII.pdf` | Hiper-raio e Efimov, 11 exercícios "Sua vez" |
| 5 | `notebooks/3_tres_corpos/03_rumo_ao_efimov.ipynb` | A torre de Efimov |
| — | `HISTORICO.md` | **O que já descobrimos.** Mais recente no topo |

---

## O código — `src/`

Agrupado por número de corpos. Cada arquivo tem cabeçalho explicando a física.

### `src/comum/` — serve para todo mundo
| Arquivo | O que faz |
|---|---|
| `constantes.py` | Constantes físicas e conversão de unidades |
| `solvers.py` | Integra a equação radial: diferença central, Numerov, e `grade()`, que alinha a borda do poço à malha |
| `incerteza.py` | Extrapolação de Richardson e os **portões de validade** (`|a|/r₀ > 10`, `kR < 0.3`, `r₀/ℓ_ho < 0.1`) |

### `src/um_corpo/` — uma partícula num potencial externo
| Arquivo | O que faz |
|---|---|
| `schrodinger.py` | Oscilador harmônico, átomo de hidrogênio, caixa |

### `src/dois_corpos/` — duas partículas que interagem (o coração do laboratório)
| Arquivo | O que faz |
|---|---|
| `potenciais.py` | Poço esférico, gaussiano, Pöschl-Teller modificado, Lennard-Jones |
| `espalhamento.py` | Comprimento de espalhamento `a` e alcance efetivo `r₀` |
| `analitico.py` | Fórmulas fechadas, para conferir o numérico |
| `ajuste.py` | Ajuste da expansão de alcance efetivo `k cot δ₀ = −1/a + ½r₀k²` |
| `feshbach.py` | Ressonância de Feshbach: comprimento de espalhamento vs campo magnético |

### `src/tres_corpos/` — física de Efimov
| Arquivo | O que faz |
|---|---|
| `trimero.py` | Equação hiper-radial, com e sem armadilha. `s₀ = 1.0062378` |
| `efimov.py` | Razões universais (22.6944 e 515.03) e o parâmetro de três corpos |

### `src/literatura/` — comparação com valores publicados
| Arquivo | O que faz |
|---|---|
| `bench_referencias.py` | Potenciais e fórmulas de artigos específicos (Aziz HFD-B, esfera mole, Jeszenszki) |
| `tabelas_artigo.py` | As tabelas do artigo que reproduzimos |

### `src/relatorios/` — geração de saída
`gerar_figuras.py`, `gerar_relatorio.py`, `gerar_resultados.py`, `gerar_inicio.py`, `diario.py`

> **Os dois jeitos de importar funcionam.** `from src import trimero` (o de sempre)
> e `from src.tres_corpos import trimero` (explícito). Nenhum notebook antigo quebrou —
> `src/__init__.py` registra os apelidos.

---

## Os dados da literatura — `referencias/`

Este é o **arquivo-mãe**: tudo que a gente compara com artigo passa por aqui.

| Onde | O que tem |
|---|---|
| `literatura.py` | **44 valores, 15 fontes, 5 potenciais, 9 convenções, 3 divergências.** Cada número com artigo, página e a incerteza de arredondamento dele |
| `BENCHMARKS.yaml` | Os mesmos valores em formato legível por máquina |
| `CONVENCOES.md` | As 9 convenções que diferem entre artigos (fator ½ no Lennard-Jones, massa reduzida vs massa, sinal de δ₀) |
| `benchmarks_csv/` | Um CSV por referência: `bench_aziz`, `bench_esfera_mole`, `bench_fase_variavel`, `bench_gauss_jeszenszki`, `bench_he_dimer`, `bench_singleto_np` |
| `fichamentos/` | Suas anotações de leitura |
| `pdfs/` | Os artigos |

**Consultar pela linha de comando:**

```bash
python referencias/literatura.py            # lista tudo
python referencias/literatura.py gauss      # busca "gauss"
python referencias/literatura.py deuteron   # ignora acento
```

**A regra:** nenhum número entra no laboratório sem fonte. Se você achar um valor
sem referência, ele é suspeito.

---

## Resultados e figuras

| Onde | O que tem |
|---|---|
| `resultados/*.csv` | Tabelas reproduzidas, ajustes, convergência, validação analítica |
| `resultados/resumo.md` | Leitura em prosa dos CSVs |
| `figuras/2_dois_corpos/` | 7 figuras: poço, potenciais na unitariedade, Lennard-Jones, soluções radiais, convergência |
| `figuras/3_tres_corpos/` | Torre de Efimov, limiar vs dimensão |
| `figuras/4_experimento_39K/` | Onde procurar Efimov no potássio-39 |

---

## Testes — `tests/`

```bash
make test           # ou:  python -m pytest tests/ -q
```

**Regra da casa:** todo teste ancora em resultado **analítico** ou **publicado**,
nunca em "o que o código deu ontem". Os três últimos de `test_trimero.py` são
testes de **regressão** — cada um documenta um bug real que produziu número
errado *sem levantar exceção*, que é o pior tipo.

---

## Trabalhar no Mac e no desktop ao mesmo tempo

Ver `COMO_TRABALHAR_NAS_DUAS_MAQUINAS.md`. Resumo:

| Máquina | Comando |
|---|---|
| Windows | dois cliques em `sincronizar.bat` |
| Mac | `./infra/sincronizar.sh` |

Rode **antes de começar** e **depois de terminar**. Sempre.

---

## `arquivo/`

Coisas aposentadas. **Nada foi apagado.** Não sobe para o GitHub (está no
`.gitignore`) — existe só nesta máquina. Ver `arquivo/LEIA-ME.md`.

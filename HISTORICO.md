# Histórico do laboratório

Registro permanente. **Nada aqui se apaga** — resultado que morreu fica
registrado com a causa da morte, porque a causa costuma valer mais que o
resultado.

Formato: mais recente primeiro.

---

## 2026-08-15 · Arrumação: tudo agrupado por número de corpos

Nenhuma física mudou. Só a arrumação — mas ela custou tempo suficiente para
merecer registro, e a decisão de compatibilidade abaixo precisa ficar escrita.

### O problema

A raiz tinha 41 entradas, com duas cópias inteiras e obsoletas do projeto
(`github_publish/`, `publicar_github/`) que eram sobra do tempo em que a
publicação no GitHub era feita copiando arquivo à mão. A cópia que estava no
GitHub tinha um `solvers.py` de 74 linhas; o da raiz tem 157 — faltava a função
`grade()`, justamente a que alinha a borda do poço à malha. Ou seja: a versão
publicada estava **numericamente pior** que a local, e ninguém tinha percebido.

### O que foi feito

* Raiz: 41 → 21 entradas. **Nada foi apagado**: tudo aposentado foi para
  `arquivo/`, com `arquivo/LEIA-ME.md` dizendo o que é cada coisa e como voltar.
* `src/` reorganizado por **número de corpos**: `comum/`, `um_corpo/`,
  `dois_corpos/`, `tres_corpos/`, `literatura/`, `relatorios/`.
* `notebooks/` e `figuras/` idem.
* `referencias/`: CSVs de benchmark isolados em `benchmarks_csv/`.
* `MAPA.md`: uma página só, onde está cada coisa.

### A decisão que importa: compatibilidade dupla

Mover módulo de lugar quebra `import`. Em vez de sair reescrevendo notebook —
inclusive os que você ainda vai escrever à mão —, `src/__init__.py` registra
apelidos em `sys.modules`. Os **dois** jeitos funcionam e vão continuar:

```python
from src.tres_corpos import trimero    # explícito, diz onde mora
from src import trimero                # curto, o de sempre
```

Verificado: **26 formas de import** que existiam no projeto, todas funcionam.
`test_lab.py` 25/25, `test_solver_unificado.py` 23/23, `test_trimero.py` 11/11
nos rápidos (os de armadilha são lentos desde sempre, não regrediram).

### Uma nota de método

O bug do `solvers.py` de 74 linhas no GitHub é o tipo de coisa que arrumação
encontra e mais nada encontra. Cópia manual de pasta é uma fonte silenciosa de
divergência de versão. Daqui em diante a única ponte entre as máquinas é o git
(`sincronizar.bat` no Windows, `infra/sincronizar.sh` no Mac), justamente porque
o git não deixa isso acontecer sem avisar.

---

## 2026-08-10 (madrugada, 2) · O RESULTADO MAIS FORTE ATÉ AGORA — duas formas fechadas

### A pergunta

Em vez de comparar energias, comparar **funções de onda**. Quão parecida é a
função de onda de um canal só com a verdadeira?

No caso não-interagente as duas são gaussianas, e a sobreposição é analítica.

```
Psi_exata     ~ exp(-1/2 [ freq_radial (xi_1x^2+xi_1y^2+xi_2x^2+xi_2y^2)
                         + freq_axial  (xi_1z^2+xi_2z^2) ])
                <- ANISOTROPICA: só aperta em z

Psi_um_canal  ~ exp(-1/2 freq_efetiva * hiper_raio^2)
                <- ISOTROPICA em 6 dimensoes, por construcao.
                   É exatamente isso que uma "dimensão efetiva" produz.
```

Definindo **fidelidade** = `|<Psi_um_canal | Psi_exata>|^2` (1 = perfeita,
0 = inútil):

### O resultado com a frequência adiabática

| razão de anisotropia | fidelidade | leitura |
|---|---|---|
| 1 | 1,0000 | perfeita (é o caso isotrópico) |
| 2 | 0,9143 | boa |
| 5 | 0,5273 | sofrível |
| 10 | 0,2325 | ruim |
| 30 | 0,0395 | **inútil** |
| 100 | 0,0042 | **inútil** |

### E se eu der a MELHOR chance possível?

Em vez da média adiabática, escolher a frequência efetiva que **maximiza** a
fidelidade — o melhor que um canal isotrópico consegue fazer:

| anisotropia | freq. adiabática | fidel. | freq. ÓTIMA | fidel. ÓTIMA |
|---|---|---|---|---|
| 10 | 5,831 | 0,2325 | 2,000 | 0,4390 |
| 100 | 57,741 | 0,0042 | 2,794 | 0,0638 |
| 1000 | 577,351 | 0,0000 | **2,976** | 0,0067 |

### AS DUAS FORMAS FECHADAS

Derivadas simbolicamente (`sympy`) e conferidas numericamente até anisotropia
100.000:

```
freq_efetiva_otima  ->  3 * freq_radial       (SATURA)
fidelidade_maxima   ->  (27/4) / anisotropia  =  6,75 / anisotropia
```

| anisotropia | freq. ótima medida | fidel. máx. medida | 6,75/anisotropia |
|---|---|---|---|
| 100 | 2,7938 | 6,376e-02 | 6,750e-02 |
| 1 000 | 2,9764 | 6,710e-03 | 6,750e-03 |
| 10 000 | 2,9976 | 6,746e-04 | 6,750e-04 |
| 100 000 | 2,9998 | 6,750e-05 | 6,750e-05 |

### O que significa

**1. A frequência efetiva ótima SATURA em 3 × freq_radial.** Por mais que se
aperte o eixo z, a melhor descrição isotrópica para de acompanhar. Ela
literalmente não consegue seguir o achatamento.

**2. A fidelidade máxima cai como 1/anisotropia.** Em anisotropia 100 — valor
experimental modesto — a **melhor descrição possível** de um canal tem 6% de
sobreposição com a função de onda verdadeira.

**Isso não é "a aproximação piora". É "a aproximação deixa de descrever o
estado".** E vale para o caso **sem interação nenhuma**, que é o mais fácil
que existe.

**O problema não é a escolha do parâmetro — é a FORMA do ansatz.** Um canal
isotrópico em 6 dimensões não representa uma função de onda achatada, e
nenhum ajuste de parâmetro conserta isso.

### Onde isso deixa a discussão

Quatro resultados independentes, todos ancorados, dizendo a mesma coisa:

| # | resultado | natureza |
|---|---|---|
| 1 | espectro confinado não é geométrico (razões diferem 63%) | numérico |
| 2 | termo desprezado satura em 71% do mantido | geométrico |
| 3 | energia de um canal erra por raiz(3) = 73,2% | **exato** |
| 4 | fidelidade máxima cai como 6,75/anisotropia | **exato** |

O nº 4 é o mais forte, porque não fala de erro percentual: fala de a descrição
**deixar de ser uma descrição**.

**Ressalva permanente:** tudo isso é o caso não-interagente. É o mais amigável
possível para a teoria de dimensão efetiva. Com interação ressonante só pode
piorar — mas isso ainda precisa ser mostrado, não afirmado.

---

## 2026-08-10 (madrugada) · O erro da aproximação de um canal, medido EXATAMENTE

### O teste, e por que ele não tem parâmetro livre

Três partículas **não-interagentes** numa armadilha anisotrópica é um problema
**exatamente solúvel** — são só osciladores harmônicos independentes. Então dá
para medir o erro da aproximação de um canal contra a verdade, sem modelo,
sem ajuste, sem nada inventado.

**A conta exata.** Nove modos (3 partículas × 3 direções): 6 no plano, 3 no
eixo comprimido.

```
estado fundamental total = 3 freq_radial + 1,5 freq_axial
centro de massa          = 1 freq_radial + 0,5 freq_axial
--------------------------------------------------------
RELATIVO (a verdade)     = 2 freq_radial + 1,0 freq_axial
```

**A conta adiabática.** Um oscilador isotrópico em 6 dimensões com a
frequência efetiva, no canal hiperangular mais baixo (momento grand-angular
zero, portanto `L = 2`):

```
E_adiabatica = 3 × freq_efetiva = 3 × sqrt((2 freq_radial² + freq_axial²)/3)
```

### O resultado

Com a frequência radial fixa em 1:

| razão de anisotropia | E exata | E adiabática | erro | razão |
|---|---|---|---|---|
| 1 | 3,000 | 3,000 | +0,0% | 1,0000 |
| 2 | 4,000 | 4,243 | +6,1% | 1,0607 |
| 5 | 7,000 | 9,000 | +28,6% | 1,2857 |
| 10 | 12,000 | 17,493 | +45,8% | 1,4577 |
| 100 | 102,000 | 173,222 | +69,8% | 1,6983 |
| 10 000 | 10 002,000 | 17 320,508 | **+73,2%** | **1,7317** |

```
limite de aperto forte:  E_adiabatica / E_exata  ->  raiz(3) = 1,732051
```

**A aproximação de UM canal superestima a energia por um fator raiz(3) —
73,2% — no limite de panqueca.** Resultado fechado, exato, sem parâmetro
livre.

E ele confirma por um caminho independente a estimativa de ontem via
flutuação da fração axial (70,7%). Dois caminhos, ~70%, e este é analítico.

**Por que isso é o argumento contra a dimensão efetiva única:** uma "dimensão
efetiva" *é* exatamente uma redução a um canal hiperangular. Este cálculo
mostra que a redução erra por raiz(3) já no caso mais simples possível — sem
interação nenhuma. Com interação ressonante só piora.

### Bug encontrado pela verificação cruzada

A verificação deu `E_0 = 2,936492` em vez de `3,000000` — erro de **−2,117%,
idêntico em todas as frequências**. Testei domínio (início e fim) e grade:
**o erro não mudava em nenhum dígito**. Logo não era numérico — era fórmula.

Refazendo a substituição de variável simbolicamente (`sympy`):

```
com hiper_raio = e^x  e  u = e^(x/2) w :

    w'' = [ freq² e^(4x) − 2E e^(2x) + L² ] w
                                        ^^^
    é + L², NÃO + (L² − 1/4)
```

O `−1/4` do termo centrífugo **cancela** com o `+1/4` que a própria
substituição `u = e^(x/2) w` produz.

**O `src/trimero.py` já estava certo** — ele usa `−s0²` diretamente, e
`L² = −s0²` para o canal de Efimov. O erro estava no script de verificação.

Depois da correção, o solver reproduz o oscilador em `d` dimensões **exato em
6 casas** nos cinco casos testados:

| dimensão | grand-angular | L | E numérico | exato |
|---|---|---|---|---|
| 6 | 0 | 2,0 | 3,000000 | 3,0 |
| 6 | 1 | 3,0 | 4,000000 | 4,0 |
| 6 | 2 | 4,0 | 5,000000 | 5,0 |
| 3 | 0 | 0,5 | 1,500000 | 1,5 |
| 3 | 1 | 1,5 | 2,500000 | 2,5 |

**Lição de método:** erro que não muda com refinamento **nunca** é numérico.
Se domínio e grade não mexem no dígito, pare de mexer no numérico e volte
para a álgebra.

---

## 2026-08-10 (noite) · A armadilha anisotrópica, e por que um canal só nunca basta

### Nomenclatura (nomes completos, sem abreviação)

| símbolo curto | nome completo | o que significa |
|---|---|---|
| `xi_1`, `xi_2` | vetores de Jacobi | as duas posições relativas, cada uma com 3 componentes (x, y, z) |
| `hiper_raio` | hiper-raio | `sqrt(xi_1^2 + xi_2^2)` — o tamanho global do triângulo de 3 átomos |
| `omega_radial` | frequência radial de aprisionamento | aperto no plano (direções x, y) |
| `omega_axial` | frequência axial de aprisionamento | aperto na direção comprimida (z) |
| `anisotropia` | razão de anisotropia | `omega_axial / omega_radial`. 1 = esférico; grande = panqueca |
| `fracao_axial` | fração axial do hiper-raio | `(xi_1z^2 + xi_2z^2) / hiper_raio^2`, entre 0 e 1 |
| `omega_efetiva` | frequência efetiva hiperradial | o que sobra depois da média adiabática |
| `parametro_tres_corpos` | parâmetro de três corpos | a parede em hiper-raio pequeno; fixa onde a torre começa |

### A derivação

Somando a energia de aprisionamento das três partículas e removendo o centro
de massa, a parte relativa é:

```
V_armadilha = (1/2) m hiper_raio^2 [ omega_radial^2 (1 - fracao_axial)
                                   + omega_axial^2      fracao_axial  ]
```

**Isotrópico** (`anisotropia = 1`): a `fracao_axial` cancela. O termo vira
puramente função do hiper-raio, e a separação hiperesférica é **exata**.

**Anisotrópico**: a `fracao_axial` **não** cancela. Ela é função dos
hiper-ângulos. **Esse é o termo que quebra a separabilidade.**

### O valor médio, verificado

Amostragem uniforme na esfera de 6 dimensões, 2.000.000 de pontos:

| grandeza | medido | previsto |
|---|---|---|
| média da `fracao_axial` | **0,333089** | 1/3 = 0,333333 |
| desvio padrão | **0,235590** | — |
| mínimo / máximo | 0,0000 / 0,9991 | — |

A média 1/3 tem razão simples: `xi_1` e `xi_2` são vetores de 3 componentes, e
num estado de momento angular total zero as três direções físicas do espaço são
equivalentes. A direção comprimida carrega, em média, um terço do tamanho.

Daí a **aproximação adiabática** (trocar a `fracao_axial` pela média):

```
omega_efetiva^2 = (2 omega_radial^2 + omega_axial^2) / 3
```

— a média quadrática sobre as três direções.

### O espectro sob aperto

Com `omega_radial` fixo em `1e-8` (comprimento no plano = 10.000 unidades de
`parametro_tres_corpos`):

| anisotropia | comprimento do oscilador | N previsto | N de Efimov ligados |
|---|---|---|---|
| 1 | 10.000 | 2,95 | 3 |
| 3 | 7.227 | 2,85 | 2 |
| 10 | 4.141 | 2,67 | 2 |
| 30 | 2.402 | 2,49 | 2 |
| 100 | 1.316 | 2,30 | 2 |

Apertar mata degraus — de cima para baixo, porque o degrau de cima é o maior e
é o primeiro a não caber.

### O RESULTADO DO DIA

A aproximação adiabática troca a `fracao_axial` pela média. Mas ela **flutua**:

```
desvio_padrao / media = 0,2356 / 0,3333 = 70,7%
```

E a razão entre o termo **jogado fora** e o termo **mantido**:

```
jogado_fora / mantido = 3 * desvio_padrao * (anisotropia^2 - 1)/(anisotropia^2 + 2)
```

| anisotropia | jogado fora / mantido |
|---|---|
| 2 | 35,3% |
| 10 | 68,6% |
| 100 | 70,7% |
| 1000 | **70,7%** |

**O termo desprezado não encolhe. Ele satura em ~71% do termo mantido, por
mais que se aperte.**

E `3 x 0,2356 = 0,707` é o limite exato — vem só da geometria da esfera de 6
dimensões, não do valor da anisotropia.

**Por que isso importa:** uma "dimensão efetiva" única *é* exatamente uma
redução a um canal só. Este cálculo mostra que essa redução carrega um erro
irredutível de ~71% no termo de acoplamento, e que ele **não melhora no aperto
forte** — que é justamente onde o experimento trabalha.

Junto com o teste de auto-similaridade (abaixo), fecha o argumento: a descrição
por um parâmetro único não pode ser exata, e agora sabemos **por quanto** e
**por quê**.

**Ressalva:** o cálculo de espectro acima é a própria aproximação adiabática.
Ele mostra a truncagem da torre, mas por construção **não** produz os
cruzamentos evitados. Para vê-los é preciso o cálculo de múltiplos canais, com
o acoplamento de 71% incluído. É o próximo passo.

---

## 2026-08-10 · Trímero no hiper-raio, e a falha da dimensão fracionária

### Construído

| arquivo | o que é |
|---|---|
| `src/trimero.py` | trímero de Efimov no hiper-raio, com e sem armadilha |
| `tests/test_trimero.py` | 11 testes de âncora, incluindo 3 de regressão |
| `src/incerteza.py` | Richardson, truncamento e **portões de validade** |
| `referencias/literatura.py` | arquivo-mãe: 44 valores com fonte e convenção |

### A transformação

Com `ρ = eˣ` e `u = e^{x/2}w`:

```
w''(x) = [ ω² e^{4x} − 2E e^{2x} − s₀² ] w(x)
```

A log-periodicidade fica manifesta; a torre, que em `ρ` vive espalhada por
décadas, em `x` vira malha uniforme. **É o que torna o problema tratável.**

E a armadilha **isotrópica** entra exata, porque `Σᵢrᵢ² = 3R_cm² + ρ²`.
(A anisotrópica — a do experimento — não entra exata. É o próximo degrau.)

### Resultados

**Torre sem armadilha** — âncora passou:

| n | `κₙ` | razão |
|---|---|---|
| 0 | 6,5376e-02 | |
| 1 | 2,8792e-03 | 22,7064 |
| 2 | 1,2687e-04 | 22,6944 |
| 3 | 5,5902e-06 | 22,6944 |

Excitados: erro **3×10⁻⁷** contra `e^{π/s₀}`. Desvio concentrado no
fundamental (+0,053%) — o estado mais fundo é o que mais sente a parede, e é o
menos universal. Se o erro estivesse espalhado, seria bug.

**Armadilha trunca a torre** — `N = (s₀/π)·ln(ℓ_ho/ρ₀)`:

| `ℓ_ho` | previsto | achado |
|---|---|---|
| 10⁴ | 2,95 | 3 |
| 10³ | 2,21 | 2 |
| 10² | 1,48 | 1 |
| 10 | 0,74 | 0 |

**A dimensão fracionária falha** — e este é o resultado do dia:

| n | `Eₙ` com armadilha | razão | `Eₙ` sem | razão |
|---|---|---|---|---|
| 0 | −2,136978e-03 | **515,58** | −2,13698e-03 | 515,58 |
| 1 | −4,144788e-06 | **839,68** | −4,14480e-06 | 515,04 |
| 2 | −4,936139e-09 | | −8,04758e-09 | |

Sem armadilha as razões são geométricas em **0,11%**. Com armadilha diferem
em **62,9%**.

**A lógica:** um `d` fixo dá um `s₀(d)` fixo, que dá uma torre geométrica —
`Eₙ/Eₙ₊₁` constante para todo `n`. O espectro confinado tem razões que
diferem por 63% entre pares consecutivos. **Nenhum `d_ef` único reproduz isso.**

**O mecanismo é transparente:** a armadilha corta em `ρ ~ ℓ_ho`. O estado
fundo (`ρ ≪ ℓ_ho`) não sente — 515,58 em ambos os casos, idêntico. O estado
raso (`ρ ~ ℓ_ho`) é espremido. **O desvio depende de `n`, e um botão único não
produz dependência em `n`.**

**Ressalva honesta:** esta é armadilha *isotrópica*, não confinamento
anisotrópico. Um defensor da dimensão fracionária diria que armadilha
isotrópica não é mudança de dimensão. O argumento ainda morde porque o
mecanismo — corte numa escala fixa, afetando estados rasos mais que fundos —
é genérico a qualquer confinamento. **O caso isotrópico é o mais amigável
possível para a teoria de `d`, e já quebra em 63%.** O anisotrópico só pode
ser pior: tem isso *mais* a mistura hiperangular.

É condição **necessária**, não suficiente. Não refuta a teoria; mostra que
ela não pode ser exata.

### Bugs encontrados (todos silenciosos)

**1 · Underflow apaga os nós.** O reescalonamento contra overflow divide o
array inteiro por `1e250`. Com armadilha forte dispara várias vezes, e a
região oscilatória do início vira **zero exato** — 57.133 zeros de 60.000
pontos. A torre sumia sem erro nenhum. *Conserto:* contar nós **durante** a
integração. Teste de regressão: `test_regressao_underflow_apaga_os_nos`.

**2 · Overflow no produto de sinais.** `w[i+1]*w[i] < 0` estoura quando ambos
são `~1e250`. *Conserto:* comparar sinais.

**3 · Colisão de nome.** `n` era o índice do estado e o tamanho da grade.
*Conserto:* grade virou `npts`.

O nº 1 é o perigoso: número errado, sem exceção, sem aviso.

---

## 2026-08-10 · Incerteza e portões de validade

`src/incerteza.py`. Richardson sobre a grade, teste de truncamento, e
**portões**: `|a|/r₀ > 10`, `kR < 0,3`, `r₀/ℓ_ho < 0,1`.

Relatório de conformidade agora sai honesto:

```
gauss deuteron a      5.4002699 ± 0.0000010 fm    ref=5.4     0.0σ  OK
poco nn a (central)  -18.5277277 ± 0.0000064 fm   ref=-18.52  1.5σ  OK
poco nn a (numerov)  -18.5277 ± 0.0063 fm         ref=-18.52  1.0σ  OK
```

**O Achado nº 2 virou barra de erro:** o Numerov no poço carrega incerteza
1000× maior, automaticamente.

**Dois defeitos que o módulo pegou em si mesmo:**

- Comparar `5.400269851907 ± 5e-11` com `5.40` dava **4.871.117σ**. O erro
  estava na comparação: quando o artigo publica `5.40`, a incerteza dele é
  `±0,005` e **domina a nossa**.
- Alertas de "ordem errada" quando o resultado já convergiu à precisão de
  máquina — ali a ordem medida é ruído.

**O portão apontado para a literatura** achou que `|a|/r₀` do dêuteron é
**3,10**. O exemplo canônico de universalidade em livro-texto não está no
regime universal. Não é erro de ninguém — fica invisível até automatizar.

E para o ³⁹K: a janela universal é `B ∈ (401,65 · 403,88) G`, **2,2 G de
largura**, porque `r_eff = 136 a₀` é 2,1× o `R_vdW`. Duas das três previsões
da feição de Efimov (404,63 e 405,27 G) caem **fora** dela.

---

## 2026-08-10 · Solver confinado de dois corpos

Âncora: **Busch, Englert, Rzażewski, Wilkens (1998)** — dois átomos numa
armadilha harmônica, alcance zero, solução analítica fechada.

Confirmado: unitariedade `E = 0,5000000000 ħω`; não-interagente `1,4999920`.

**Resultado:** na unitariedade, `ΔE = 0,2957 (r₀/ℓ_ho) ħω` no limite
`r₀ → 0`. Linear e limpo sobre um fator 8 em `r₀`.

Minha receita ingênua de comprimento de espalhamento dependente de energia
previa `0,3989` — **erra 26%**. A prescrição correta está em Blume & Greene
(2002) e Bolda–Tiesinga–Julienne (2002). **Não chutar de novo.**

**Erro de fundo que eu cometi:** comparei com Busch em pontos onde `|a|/r₀`
valia 0,6 a 4. Busch é alcance zero e exige `|a| ≫ r₀`. A "discordância" que
parecia física era hipótese violada. **Foi esse erro que gerou o
`incerteza.py`.**

**Aplicação:** a fronteira onde alcance finito passa a importar no ³⁹K fica em
`ℓ_z ≈ 72 nm`, ou `ν_z ≈ 100 kHz`. *Pergunta para a Patrícia: qual `ν_z` o
acordeão alcança no aperto máximo?*

---

## 2026-08-09 · `ρ_d` — resultado que morreu, e o que sobrou

Tentativa: generalizar a integral de Bethe (`r₀`) para `d` dimensões, para
ver como o alcance efetivo se comporta no crossover.

Âncoras em `d=3` passaram: poço `ρ₃/R = 0,99983` (exato `R`), mPT
`ρ₃·μ = 2,000000` (exato `2/μ`).

**Morreu em dois testes:**

1. O `ρ` do mPT depende de onde eu corto a cauda — de `−0,47` (R=8) a `−9,16`
   (R=24). Física não depende de truncamento. A integral **não converge** para
   `d < 3`.
2. `ρ/R` não é invariante de escala fora de `d=3`.

**O que sobrou, e vale:** `ρ_d ∝ R^{4−d}`, exato em **cinco casas** para todo
`d` testado. Ou seja, **`ρ_d` não tem dimensão de comprimento exceto em
`d=3`.**

Logo a pergunta "como `r₀` varia no crossover" está **mal-posta**. Não existe
`r₀(d)` contínua ligando 3D a 2D, porque em 2D a expansão de baixa energia é
logarítmica e não tem termo `−1/a`. **Muda a forma da parametrização, não só
o valor dos parâmetros.**

Isso alimenta a hipótese testada em 2026-08-10 e confirmada.

---

## Pendências

| # | o que | quem |
|---|---|---|
| 1 | `.git/config.lock` — apagar à mão | Pedro |
| 2 | `github_publish/` obsoleto — arquivar | Pedro |
| 3 | DOI do Zenodo | Pedro |
| 4 | `ν_z` máximo do acordeão, densidade, `L₃` | perguntar à Patrícia |
| 5 | Equação transcendental em `d` fracionário (arXiv:1708.00012) | ler |
| 6 | Prescrição correta de alcance finito em armadilha (Blume & Greene 2002) | ler |
| 7 | Armadilha **anisotrópica** no hiper-raio | próximo passo de código |
| 8 | Ponte para `L₃` (Braaten–Hammer) | falta |
| 9 | QMC (VMC/DMC) | não começou |

## Leitura, em ordem

1. **Etrych et al., PRR 5, 013174 (2023)** — aberto, curto, sistema da Patrícia
2. **Madeira et al., PRA 104, 033301 (2021)** — o método da dissertação
3. **Levinsen, Massignan, Parish, PRX 4, 031020 (2014)** — quase-2D de verdade
4. *Nonuniversal EoS of a Quasi-2D Bose Gas*, arXiv:2402.04703 (2024)
5. Madeira, Few-Body Syst. 65, 70 (2024) — o fio do alcance efetivo

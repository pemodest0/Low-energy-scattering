# Visita ao laboratório da Patrícia — ³⁹K

Pedro Modesto · IFSC-USP · agosto de 2026

> **Como usar:** não é apresentação. É uma lista de números que eu preciso trazer de
> volta, porque **cada um deles muda o que eu calculo**. Se eu sair de lá com a
> tabela da Seção 2 preenchida, a visita valeu.

---

## 1. O que eu levo

Uma figura e três números. Só isso — primeira visita não é seminário.

**`figuras/39K_onde_procurar_efimov.png`**

A curva `a(B)/R_vdW` para a ressonância larga do |1,1⟩, com **três previsões
diferentes** de onde deve estar a primeira feição de perda de Efimov:

| origem | `a₋` (a₀) | campo previsto |
|---|---|---|
| Zaccanti 2009 (medido) | −1500 | **403,77 G** |
| Etrych 2023 (medido) | −830 ± 40 | **404,63 G** |
| universalidade vdW (`−9,7 R_vdW`) | −627 | **405,27 G** |

Elas estão **1,5 G separadas**. Isso é o argumento que eu levo: a divergência da
literatura vira uma diferença de campo que um experimento consegue resolver.

Os parâmetros que usei: `B₀ = 402,74 G`, `Δ = −52,1 G`, `a_bg = −29 a₀`,
`R_vdW = 64,6 a₀`. Confirmar se batem com os dela.

---

## 2. O que eu preciso perguntar

Ordenado por quanto muda o meu cálculo. As três primeiras são as que importam.

### 2.1 O confinamento — define o eixo inteiro da minha dissertação

| pergunta | por que muda o meu cálculo | resposta |
|---|---|---|
| Qual a faixa de espaçamento da rede em acordeão? (mínimo e máximo) | define o `ℓ_z` acessível, e portanto **qual pedaço do crossover 3D→2D existe no laboratório** | |
| `ω_z` correspondente, e `ω_r` no plano? | `ℓ_z = √(ħ/mω_z)` é a escala que entra na minha conta | |
| A rede muda de espaçamento **durante** um experimento, ou entre experimentos? | crossover dinâmico ou quase-estático — são físicas diferentes | |
| Qual `ℓ_z/a` vocês conseguem alcançar? | é o parâmetro adimensional do crossover. Preciso saber se dá para chegar em `ℓ_z ≲ a` | |

### 2.2 A ressonância — define onde eu calculo

| pergunta | por que | resposta |
|---|---|---|
| Vocês usam a larga do \|1,1⟩ em ~402,7 G, ou alguma estreita? | a larga tem `s_res = 2,8` — **intermediária**, não larga de verdade. Isso quebra a universalidade vdW (é a nossa divergência D2) | |
| Qual `B₀`, `Δ` e `a_bg` vocês usam na prática? | para eu corrigir o meu `feshbach.py` | |
| Qual a estabilidade e a resolução do campo? | define se 1,5 G é distinguível — e portanto se a pergunta da Seção 1 é respondível | |

### 2.3 A medida — define o observável que eu tenho que prever

| pergunta | por que | resposta |
|---|---|---|
| Vocês medem `L₃` (taxa de perda de três corpos)? Com que precisão? | é **o** observável de Efimov. Se medem, eu tenho alvo | |
| Densidade e temperatura típicas | `L₃` depende de `n²`; e a temperatura corta a feição (efeito térmico borra o pico) | |
| Número de átomos e estabilidade tiro a tiro | define o contraste mínimo detectável numa feição de perda | |
| Tempo de vida do gás no campo de interesse | perto da ressonância a nuvem morre rápido — quanto tempo há para medir? | |
| Medida in situ ou depois de tempo de voo? | muda o que eu comparo | |

### 2.4 Cronograma

- Quando a rede em acordeão entra no experimento? (o projeto diz "até o fim de 2026")
- O que já está pronto e o que falta?
- Existe alguma medida **já feita** que eu possa tentar reproduzir agora, como
  calibração do meu código? (mesmo que seja 3D puro)

---

## 3. A pergunta forte

Se der abertura, essa é a que vale a visita inteira:

> Zaccanti (2009) e Etrych (2023) medem `a₋` na **mesma ressonância do mesmo átomo**
> e discordam por um fator ~1,8 (−1500 contra −830 a₀). Isso é conhecido no grupo?
> O experimento de vocês teria condições de arbitrar?

É uma divergência aberta na literatura, no sistema exato dela. Se a resposta for
"dá para medir", isso é um resultado publicável e minha dissertação tem um alvo
experimental concreto desde o primeiro ano.

---

## 4. O que eu posso oferecer

Sem prometer o que não tenho.

**Já funciona hoje:**

- `a(B)` e a posição prevista das feições, com os parâmetros que ela escolher
- `a` e `r₀` de potenciais microscópicos, validados contra literatura em 4 camadas
  (25 testes, erro 10⁻⁵–10⁻¹³)
- o `r₀` — que é o meu diferencial: quase toda a teoria de Efimov publicada é de
  **alcance zero**, e o `r_eff = 136 a₀` que o Etrych mede não é pequeno comparado
  ao `R_vdW = 64,6 a₀`

**Em construção (não prometer prazo):**

- o espectro de trímeros sob confinamento, via Monte Carlo quântico, estendendo
  Madeira *et al.* PRA **104**, 033301 (2021)
- previsão de `L₃(B, ℓ_z)` — o que ela mediria ao girar a rede

---

## 5. Depois da visita — o que fazer com as respostas

| se ela disser | eu faço |
|---|---|
| a faixa de `ℓ_z` é X | rodo `s₀(d)` só na janela de `d` que o experimento alcança |
| usam a ressonância de 402,7 G | atualizo `feshbach.py` para `B₀ = 402,74` e uso `r_eff = 136 a₀` |
| medem `L₃` com precisão Y | calculo o contraste esperado da feição e digo se é detectável |
| a rede atrasou | uso a folga para fechar o QMC 3D primeiro, reproduzindo o PRA 2021 |

---

## Referências que eu devo ter no bolso

- Etrych *et al.*, **PRR 5, 013174 (2023)** — os números modernos do ³⁹K. Aberto.
- Zaccanti *et al.*, **Nature Physics 5, 586 (2009)** — o espectro original.
- Levinsen, Massignan, Parish, **PRX 4, 031020 (2014)** — trímeros de Efimov sob
  confinamento forte. É o cálculo de quase-2D de verdade.
- Sandoval *et al.*, **J. Phys. B 51, 065004 (2018)**, arXiv:1708.00012 —
  "Squeezing the Efimov effect". Dimensão fracionária; Yamashita e Frederico são
  coautores do Lucas no PRA 2021.
- Madeira, **Few-Body Syst. 65, 70 (2024)** — o papel do alcance efetivo.

Tudo com fonte, convenção e condição de validade em `referencias/literatura.py`.

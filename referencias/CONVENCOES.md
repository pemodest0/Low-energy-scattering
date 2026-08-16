# Convenções — a pedra de Roseta

> **Regra do repositório:** nenhum valor entra em `BENCHMARKS.yaml` sem
> apontar para um `id` deste arquivo. Número sem convenção declarada é
> número sem significado.
>
> Motivo: os dois erros mais caros deste projeto — o fator 2 do
> Lennard-Jones e o ℓ_vdW que saiu 54,3 em vez de 64,6 a₀ — **não** foram
> falta de referência. Foram convenção não declarada.

---

## C1 · Prefator do potencial e qual massa

Na equação radial só entra a combinação

```
U(r) ≡ 2 m_r V(r) / ħ²        u'' = [U − k² + l(l+1)/r²] u
```

Três escritas circulam na literatura:

| id | escrita de V | U resultante |
|---|---|---|
| `C1a` | V = (ħ²/2m_r)·W | U = W |
| `C1b` | V = (ħ²/m_r)·W | U = 2W |
| `C1c` | V = (ħ²/m)·W (massa **cheia**) | U = W, se m = 2m_r |

⚠️ Para partículas **idênticas** m_r = m/2, logo `C1b` e `C1c` diferem por
exatamente 2. É a origem da maior parte dos erros de fator 2 do campo.

**Neste repositório:** o núcleo é adimensional com ħ = m_r = 1, logo
U = 2V — ou seja, escrevemos os potenciais na convenção **`C1b`**.

**No artigo RBEF 45, e20230079:** convenção `C1b` para poço, mPT e
gaussiano (verificado). A Eq. (121), do Lennard-Jones, está impressa em
`C1b` mas as constantes da Tabela 4 **só reproduzem em `C1a`** — ver
`DIVERGENCIAS` abaixo.

## C2 · Escala de van der Waals

Três comprimentos diferentes, frequentemente com o mesmo nome:

| id | símbolo | definição | ³⁹K |
|---|---|---|---|
| `C2a` | R_vdW | ½·(2μC₆/ħ²)^{1/4} (Chin et al., RMP 82, 1225) | 64,6 a₀ |
| `C2b` | ā | comprimento **médio** de espalhamento, ā = [4π/Γ(¼)²]·R_vdW ≈ 0,9560·R_vdW | ~61,8 a₀ |
| `C2c` | — | definições sem o fator 2 na massa | 54,3 a₀ ← **erro já cometido** |

⚠️ O `μ` de `C2a` é a **massa reduzida**, e o fator 2 dentro da raiz é
parte da definição, não da derivação. Conferir sempre a definição, não a
fórmula.

## C3 · Sinal do comprimento de espalhamento

- **Convenção atômica (a nossa):** ψ ~ 1 − a/r; a > 0 significa
  repulsivo-equivalente / estado ligado presente.
- **Convenção nuclear de nêutrons:** o comprimento de espalhamento
  coerente tabelado (`b`) é frequentemente definido com o sinal oposto,
  a partir de ψ ~ 1 − b/r **com a onda espalhada escrita como −b·e^{ikr}/r**.

⚠️ Tabelas de `b` para nêutrons térmicos (as que se usa em difração e em
física médica) **podem** trazer o sinal invertido em relação a este
repositório. Ao importar qualquer valor de tabela de nêutrons, conferir o
sinal contra um caso conhecido (o n-p tripleto deve dar a > 0).

## C4 · Expansão de alcance efetivo

```
k·cot δ₀ = −1/a + ½ r₀ k² + O(k⁴)        ← a nossa (C4a)
```
Variantes na literatura trocam o sinal de r₀ ou absorvem o ½. Ao importar
r₀, verificar que o limite de alcance zero dá r₀ → 0 e que o dêuteron dá
r₀ > 0.

## C5 · Unidades de energia

Nunca implícitas. Fatores usados aqui (CODATA 2018, `src/constantes.py`):

| grandeza | valor |
|---|---|
| ħc | 197,3269804 MeV·fm |
| k_B | 1,380649×10⁻²³ J/K |
| a₀ | 0,529177210903 Å |
| 1 Å | 10⁵ fm |

Energias de gases frios aparecem como **mK, µK, nK ou h·Hz**. Conversão
explícita sempre.

## C6 · Parâmetro de três corpos

A mesma física em quatro sotaques:

| id | símbolo | o que é |
|---|---|---|
| `C6a` | a₋ | valor de a (negativo) onde o trímero toca o limiar de 3 átomos livres |
| `C6b` | κ* | número de onda do trímero na unitariedade |
| `C6c` | R₀ | corte de curto alcance no hiper-raio (o que o nosso `efimov.py` usa) |
| `C6d` | Λ | cutoff em espaço de momento (EFT) |

⚠️ a₋ tem **índice de nível**: a₋⁽¹⁾ é o primeiro trímero, a₋⁽²⁾ o segundo
(≈ 22,7× maior). Citar sem o índice gera confusão — e valores publicados
de "a₋" para o mesmo átomo podem se referir a níveis diferentes ou a
**tetrâmeros**. (Erro já cometido: −690 a₀ citado de memória para o ³⁹K,
quando −1500 a₀ é o trímero e ≈ −650 a₀ é o tetrâmero associado.)

## C7 · Razões de Efimov

| grandeza | valor | observação |
|---|---|---|
| s₀ | 1,0062378 | três bósons idênticos, alcance zero |
| e^{π/s₀} | 22,694 | razão de **comprimento** (e de κ) |
| e^{2π/s₀} | 515,03 | razão de **energia** |

⚠️ "fator de Efimov" pode significar 22,7 **ou** 515 conforme o autor.
Sempre dizer de qual grandeza.

## C8 · Alcance do potencial vs alcance efetivo

`R` (onde |V| < ε, uma escolha numérica) **não é** `r₀` (a integral da
Eq. 56). Coincidem apenas na unitariedade, para o poço esférico. Não
misturar ao comparar com literatura.

---

# DIVERGÊNCIAS REGISTRADAS

> Divergência se registra, não se apaga.

### D1 · Fator 2 na Eq. (121) do RBEF 2023 (Lennard-Jones)

Convenção declarada: `C1b`. Convenção que reproduz a Tabela 4: `C1a`.

Teste de controle, caso nêutron-nêutron, mesma convenção aplicada aos três:

| potencial | em `C1b` | em `C1a` | publicado |
|---|---|---|---|
| poço esférico | **−18,528** ✓ | — | −18,52 |
| gaussiano | **−18,554** ✓ | −1,685 | −18,55 |
| Lennard-Jones | +2,055 ✗ | **−18,499** ✓ | −18,5 |

Dois potenciais obedecem, um não → descarta erro de leitura global.
Efeito: dobrar C₁₂ e C₆ juntos mantém a posição do mínimo
(r_m⁶ = 2C₁₂/C₆) e **dobra a profundidade** (U_min = −C₆²/4C₁₂); em `C1b`
o n-n adquire um nó, i.e. o dinêutron existiria.

Correções equivalentes: (i) Eq. 121 com ħ²/2m_r, ou (ii) constantes da
Tabela 4 divididas por 2.

**Status:** a levar ao Lucas (coautor). Ver `PAUTA_LUCAS.md`.
**Nosso código:** implementa `C1a` para o LJ (ver nota na classe
`LennardJones`), portanto reproduz a Tabela 4.

### D2 · a₋ do ³⁹K vs universalidade de van der Waals — **RESOLVIDA**

**Status: não era erro nosso. É física publicada.**

Nossa previsão universal estava **certa**: Etrych et al. 2023 escrevem,
com todas as letras, que a universalidade de Efimov–van der Waals daria
`a₋ = −9,7·r_vdW ≈ −630 a₀` para o ³⁹K. É exatamente o número que a
nossa conta dava.

O que a natureza faz é **outra coisa**. Medindo em quatro ressonâncias:

| estado | B_res (G) | s_res | r_eff (a₀) | a₋ (a₀) | a₋/r_vdW | η₋ |
|---|---|---|---|---|---|---|
| \|1,1⟩ | 402,74(1) | 2,8 | 136 | **−830(40)** | −12,8(6) | 0,27(6) |
| \|1,0⟩ | 472,33(1) | 2,8 | 137 | −840(30) | −13,0(5) | 0,26(5) |
| \|1,−1⟩ | 33,5820(14) | 2,6 | 135 | −908(11) | −14,05(17) | 0,25(1) |
| \|1,−1⟩ | 162,36(2) | 1,1 | 59 | −780(70) | −12,1(11) | 0,5(1) |
| \|1,−1⟩ | 561,14(1) | 2,5 | 132 | −810(30) | −12,5(5) | 0,33(6) |

Conclusão dos autores: **quebra ubíqua da universalidade de Efimov–van
der Waals no ³⁹K**, com a₋ = −13(1)·r_vdW em vez de −9,7·r_vdW.

Mas — e isto é o ponto fino — as **razões** de Efimov (posições
relativas: a_p/a₋ = −1 universal; medido −1,08(9)) **continuam
obedecendo** a teoria. Ou seja: a torre tem a escala certa, mas está
**deslocada**. O que falha é o parâmetro de três corpos, não a
invariância de escala discreta.

Origem provável: estas são ressonâncias de **força intermediária**
(s_res ≈ 1–3). A universalidade vdW foi entendida para ressonâncias
largas (s_res ≫ 1), via reflexão quântica no potencial de van der Waals;
para estreitas, quem fixa a₋ é a força da ressonância. No meio, "the
situation is more open" — palavras dos autores.

✅ **Registro do diário CONFIRMADO na fonte primária.** O valor
`a₋ = −1500(90) a₀` está corretamente citado: Zaccanti et al. 2009
(arXiv:0904.4453) escreve textualmente *"A fit of the main maximum with
equation (2) gives a₁ = −1500(90) a₀, and η = 0.14(2)"*, para a
ressonância larga do ³⁹K no estado \|1,1⟩ — **a mesma** que o Etrych
mede. Ver DIVERGÊNCIA D3.

**Consequência prática:** a previsão do 1º trímero do ³⁹K muda. Com
a₋ = −830 a₀ em vez de −1500 a₀, o campo magnético do primeiro pico de
perda se desloca. Recalcular com o `feshbach.py`.

---

## C9 · Energia do dímero perto da ressonância (forma refinada)

A forma ingênua `E_b = ħ²/(2m_r a²)` tem uma correção conhecida de
alcance. Etrych 2023, Eq. (2), usa

```
E_b = ħ² / [ m (a − ā)² ] ,     ā = 0,956·r_vdW
```

com **m = massa do átomo** (não a reduzida). Note que ħ²/(m a²) e
ħ²/(2m_r a²) são a mesma coisa para partículas idênticas — a diferença
real é a **subtração de ā**, que desloca o polo. Para o ³⁹K,
ā ≈ 0,956 × 64,6 ≈ 61,8 a₀.

Usar esta forma ao comparar E_b(B) com dados experimentais.


### D3 · Zaccanti 2009 vs Etrych 2023 — dois experimentos, uma ressonância

**Status: discrepância aberta na literatura. Não é erro de registro nosso.**

Mesma espécie, mesmo estado (\|1,1⟩), mesma ressonância larga (~402,7 G),
mesmo observável (máximo de K₃ em a < 0):

| fonte | a₋ (a₀) | a₋/r_vdW | η₋ | método |
|---|---|---|---|---|
| Zaccanti 2009 | **−1500(90)** | −23,2 | 0,14(2) | armadilha harmônica, K₃(a) |
| Etrych 2023 | **−830(40)** | −12,8(6) | 0,27(6) | caixa homogênea, T < 100 nK, espectroscopia de estado ligado |

Fator ≈ 1,8 entre os dois.

**Hipótese testada e DESCARTADA — calibração de campo.** As duas fontes
usam B₀ diferentes (402,50 vs 402,74 G). Testamos numericamente
(`src/feshbach.py`):

- No campo onde o Etrych vê −830 a₀ (B = 404,650 G), a calibração do
  Zaccanti daria **−732 a₀** — precisaria de 2,05× para chegar a −1500.
- No campo onde o Zaccanti vê −1500 a₀ (B = 403,527 G), a calibração do
  Etrych daria **−1973 a₀** — a diferença *piora*.

A discrepância não é de conversão B → a.

**Pistas dentro do próprio Zaccanti 2009** (os autores as levantam):
1. η₋ = 0,14(2) é *"unexpectedly larger"* que o valor medido do lado
   a > 0 — e η deveria ser constante no espectro universal;
2. eles suspeitam de uma **ressonância de quatro corpos em 0,90·a₁**
   alargando o pico de três corpos;
3. *"a close inspection of the region around −1500 a₀ reveals a doublet
   structure"* — dois picos parcialmente resolvidos.
4. Eles veem **dois** máximos: −1500 e −650 a₀, e atribuem o de −1500 ao
   trímero e o de −650 ao tetrâmero associado.

**Observação nossa (não é afirmação):** o valor do Etrych, −830 a₀, cai
*entre* os dois picos do Zaccanti. E o pico de −650 a₀ do Zaccanti é
notavelmente próximo da previsão universal de −630 a₀. Se a atribuição
trímero/tetrâmero do Zaccanti estivesse trocada ou contaminada, os dois
experimentos se aproximariam. **Não temos como decidir isso daqui.**

**Argumento a favor do Etrych:** consistência interna. Ele mede
a₋/r_vdW = −13(1) em **quatro** ressonâncias independentes, com
espectroscopia de estado ligado para fixar B_res em 10 mG.

**Pergunta para o Lucas e para a Patrícia.** Se a resposta não for óbvia
para eles, isto é material de discussão real.

**Decisão de código:** `src/feshbach.py` mantém os dois valores, com a
fonte declarada. Previsão do 1º trímero: **403,53 G** (Zaccanti) ou
**404,65 G** (Etrych) — ~1,1 G de diferença, resolvível no experimento.

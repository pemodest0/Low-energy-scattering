# Etrych et al. 2023 — Pinpointing Feshbach resonances and testing Efimov universalities in ³⁹K

**Referência:** Phys. Rev. Research **5**, 013174 (2023)
**DOI:** 10.1103/PhysRevResearch.5.013174 · **Acesso:** aberto (CC-BY)
**Grupo:** Cavendish (Hadzibabic, Eigen) + Durham (Hutson)

## Em uma frase
Mapeia com alta precisão 8 ressonâncias de Feshbach intraestado do ³⁹K e
mede as assinaturas de Efimov em quatro delas — encontrando **quebra
ubíqua da universalidade de Efimov–van der Waals**.

## Convenções que ele usa
| item | convenção | id |
|---|---|---|
| a(B) | a_bg[1 − Δ/(B − B_res)] | C5 |
| E_b do dímero | ħ²/[m(a − ā)²], **m = massa do átomo** | **C9** |
| ā | 0,956·r_vdW; r_vdW = 64,6 a₀ para ³⁹K | C2a / C2b |
| a₋ | valor de a onde o trímero **mais baixo** encontra o contínuo | C6a |

## Números extraídos → BENCHMARKS.yaml
| grandeza | valor | onde | id |
|---|---|---|---|
| B_res \|1,1⟩ | 402,74(1) G | Tab. I | `k39_res_larga_Bres` |
| a_bg·Δ | 1530(20) a₀·G | Tab. I | `k39_res_larga_abgDelta` |
| B_zero | 350,4(1) G | Tab. I | `k39_res_larga_Bzero` |
| s_res | 2,8 | Tab. III | `k39_s_res_402` |
| r_eff | 136 a₀ | Tab. III | `k39_reff_402` |
| a₋ | −830(40) a₀ | Tab. III | `k39_a_menos_402` |
| a₋/r_vdW | −12,8(6) | Tab. III | — |
| a_p/a₋ | −1,08(9) | Tab. III | `k39_razao_ap_amenos` |
| a₋ universal | −9,7 r_vdW ≈ −630 a₀ | legenda Tab. III | `k39_a_menos_universal` |

## O achado central
**O parâmetro de três corpos quebra a universalidade; as razões de Efimov
não.**

- a₋/r_vdW medido: −13(1) em quatro ressonâncias (vs −9,7 universal)
- a_p/a₋ medido: −1,08(9) (vs −1 universal) → **obedece**

Leitura: a torre de Efimov tem a **escala certa** (invariância de escala
discreta intacta) mas está **deslocada** (âncora de curto alcance
diferente do previsto). Estas são ressonâncias de força **intermediária**
(s_res ≈ 1–3); a universalidade vdW foi entendida para largas (s_res ≫ 1),
via reflexão quântica no potencial de vdW.

## O que disso vira código
1. Atualizar `RESSONANCIAS_39K` com B_res = 402,74(1) e Δ de a_bg·Δ = 1530.
2. Implementar C9 (E_b com a − ā) no `feshbach.py` — o E_b(B) fica
   comparável direto com o experimento.
3. Refazer a previsão do 1º trímero com a₋ = −830 a₀ (ver abaixo).

## Consequência imediata para o nosso laboratório
Previsão do 1º pico de perda de Efimov no ³⁹K, ressonância larga:

| a₋ usado | origem | B previsto |
|---|---|---|
| −1500 a₀ | diário (a reconferir) | 403,53 G |
| **−830 a₀** | **Etrych 2023** | **404,39 G** |
| −630 a₀ | universalidade vdW | 405,01 G |

(com B₀ = 402,74 e Δ = −52,76 de Etrych: **404,65 G**)

## O que eu não entendi
- Por que ressonância de força intermediária desloca a₋ para MAIS
  negativo, e não para menos? Qual é o mecanismo físico?
- Como η₋ (parâmetro de largura de Efimov) se relaciona com s_res?

# -*- coding: utf-8 -*-
"""
Gera INICIO.html — hub único de estudo: todos os guias, resultados e
figuras embutidos num só arquivo (sem dependências, abre no navegador).

Rode:  python -m src.gerar_inicio
"""
import base64
import os

import markdown

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _b64(caminho):
    with open(caminho, "rb") as f:
        return base64.b64encode(f.read()).decode()


def _img(rel, legenda, olhar):
    caminho = os.path.join(AQUI, rel)
    if not os.path.exists(caminho):
        return f"<p><em>(figura {rel} não encontrada — rode main.py)</em></p>"
    return (f'<figure><img src="data:image/png;base64,{_b64(caminho)}" '
            f'alt="{legenda}">'
            f'<figcaption><b>{legenda}</b><br>👁️ <em>O que olhar: {olhar}</em>'
            f'</figcaption></figure>')


def _md(rel):
    caminho = os.path.join(AQUI, rel)
    if not os.path.exists(caminho):
        return f"<p><em>({rel} não encontrado)</em></p>"
    texto = open(caminho, encoding="utf-8").read()
    return markdown.markdown(texto, extensions=["tables", "fenced_code"])


def comece_aqui():
    passos = """
<h2>🚀 Só 3 coisas para saber</h2>
<ol class="passos">
<li><b>Estudar lendo (agora, aqui mesmo):</b> role esta aba — a física e
os gráficos estão logo abaixo, com "o que olhar" em cada um. As outras
abas têm os guias completos. <b>Você não precisa abrir mais nada.</b></li>
<li><b>Estudar interagindo:</b> dois cliques em <code>app.bat</code> abrem
o AMBIENTE UNIFICADO (7 estações: simular + teoria + código + resultados).
Para os notebooks, <code>estudar.bat</code> — comece pelo
<code>00_do_zero</code> (do NADA, devagar, uma ideia por célula); depois o
<code>01_visita_guiada</code>; depois o <code>02_formulas_na_pratica</code>
(cada fórmula com tradução símbolo a símbolo e um "🔧 mexa aqui") e o
<code>03_rumo_ao_efimov</code> (a ponte para o seu mestrado: torre de
Efimov e o ³⁹K de verdade). Ou <code>laboratorio.bat</code> para os
sliders.</li>
<li><b>Quando quiser profundidade:</b> abas "Guia" (física + pegadinhas) e
"Repositório" (o que estudar a fundo). O podcast do NotebookLM combina
com esta página.</li>
</ol>

<h2>A física em 3 frases</h2>
<p>Em baixa energia, tudo que um potencial faz se resume a dois números:
o <b>comprimento de espalhamento a</b> e o <b>alcance efetivo r₀</b>.
Fora do alcance do potencial, a função de onda de energia zero é sempre a
reta <code>g₀(r) = 1 − r/a</code>; o intercepto com o eixo é o próprio a.
Sinal de a conta a história: a&lt;0 = quase liga (nêutron-nêutron),
|a|→∞ = limiar exato (unitariedade, o regime do efeito Efimov do seu
mestrado), a&gt;0 = estado ligado raso (dêuteron, 1 nó).</p>


<h2>📐 As fórmulas do laboratório (tradução em palavras)</h2>
<p>Cada uma destas vive numa célula do notebook <code>02_formulas_na_pratica</code>
com código pra rodar. Aqui vai a versão de bolso:</p>

<div class="formula"><code>u''(r) = 2·V̄(r)·u(r)</code>
<p><b>A equação-mãe</b> (radial, energia zero, onda-s). Lê-se: <em>a curvatura
de u é o potencial vezes u</em>. Onde V̄&lt;0 a curva entorta pra baixo;
onde V̄=0 a curvatura zera e u vira uma RETA. u(r) = r·ψ(r) é a função de
onda "reduzida"; r é a distância ENTRE as duas partículas (fm).</p></div>

<div class="formula"><code>g₀(r) = 1 − r/a</code>
<p><b>A reta universal.</b> É a cara de u fora do alcance. Cruza zero em
r = a: o comprimento de espalhamento é literalmente o intercepto de uma
reta. a&lt;0 → cruza "atrás da origem" (não liga); a→∞ → reta horizontal
(unitariedade); a&gt;0 → cruza na frente (estado ligado, 1 nó).</p></div>

<div class="formula"><code>a = R − 2Δr·u(R) / [u(R+Δr) − u(R−Δr)]</code>
<p><b>Como o computador lê o a (Eq. 110).</b> Pura geometria de reta:
valor u(R) dividido pela inclinação (estimada pelos dois vizinhos, a
"diferença central") diz onde a reta cruza zero. Δr é o passo da grade.</p></div>

<div class="formula"><code>r₀ = 2·∫₀ᴿ [g₀²(r) − u₀²(r)] dr</code>
<p><b>O alcance efetivo (Eq. 56).</b> Área entre o "mundo ideal" (a reta g₀,
potencial de alcance zero) e o "mundo real" (u₀) — mede o tamanho da região
onde o potencial de fato age. Antes, u é normalizada para colar na reta em
R: C = (1 − R/a)/u(R) (Eq. 111).</p></div>

<div class="formula"><code>a = R·[1 − tan(√2v̄)/√2v̄]</code>
<p><b>O poço esférico exato (Eq. 80).</b> A tangente manda: quando
√2v̄ = π/2 (v̄ = π²/8 ≈ 1,2337) ela explode e a→±∞ — nasce o primeiro
estado ligado. Cada nova explosão (3π/2, 5π/2…) = mais um estado. É a
fórmula por trás das "paredes" da Fig. 10.</p></div>

<div class="formula"><code>a·μ = (π/2)·cot(πλ/2) + γ + Ψ(λ),&nbsp; v̄ = λ(λ−1)/2</code>
<p><b>O Pöschl-Teller exato (Eq. 117).</b> γ ≈ 0,5772 é a constante de
Euler-Mascheroni e Ψ é a função digamma. Na unitariedade λ=2 (v̄=1) e
r₀ = 2/μ exato — por isso o mPT é o potencial favorito dos QMCs.</p></div>

<div class="formula"><code>E_zr = −ħ²/(2m_r·a²) &nbsp;&nbsp;|&nbsp;&nbsp; 1/a = κ − r₀κ²/2 → E_fr = −ħ²κ²/(2m_r)</code>
<p><b>De dois números à energia (Eqs. 93/95).</b> Só com a: E_zr = −1,416 MeV
(64% do dêuteron real). Somando r₀: E_fr = −2,223 MeV vs −2,224 MeV medido.
κ é a taxa com que a onda ligada morre: ψ ~ e^(−κr).</p></div>

<div class="formula"><code>u_{i+1} = 2u_i − u_{i−1} + 2(Δr)²·V̄_i·u_i</code>
<p><b>Diferença central (Eq. 99)</b> — a marcha de dominós: dois pontos
anteriores dão o próximo. Erro ~(Δr)². O <b>Numerov (Eq. 101)</b> adiciona
correções [1 ± h²ξ/12] que cancelam erros até (Δr)⁴ — MAS pressupõe
potencial suave: no poço com borda abrupta ele perde do método simples
(veja o gráfico de convergência abaixo).</p></div>

<h2>Os gráficos, um a um — com cada eixo explicado</h2>"""
    figs = [
        ("figuras/fig7_potenciais_unitariedade.png",
         "Fig. 7 — os três potenciais atrativos na unitariedade",
         "EIXO X: distância r entre as duas partículas, em unidades de r₀ "
         "(r/r₀ = 1 é 'um alcance efetivo de distância'). EIXO Y: o potencial "
         "V̄ adimensional; 0 = sem interação, negativo = atração. CURVAS: azul "
         "= poço (parede abrupta em r=R — é essa quina que estraga o Numerov); "
         "laranja = Pöschl-Teller (1/cosh², decai suave); verde = gaussiana. "
         "As três têm PROFUNDIDADES e LARGURAS diferentes, mas todas foram "
         "sintonizadas para |a|=∞ e r₀=1 — mesma física de baixa energia. "
         "SE MUDAR: potencial mais fundo → sai da unitariedade pro lado a>0; "
         "mais raso → lado a<0."),
        ("figuras/fig8_lennard_jones.png",
         "Fig. 8 — Lennard-Jones na unitariedade",
         "Mesmos eixos da Fig. 7. A curva sobe a +∞ perto de r=0: é o CAROÇO "
         "REPULSIVO (C₁₂/r¹², os elétrons não deixam os átomos se sobrepor) e "
         "desce num vale raso (−C₆/r⁶, atração de van der Waals) antes de "
         "morrer. Mesmo com essa forma nada-a-ver com as da Fig. 7, os mesmos "
         "(a, r₀) saem — universalidade. SE MUDAR: ↑C₆ aprofunda o vale "
         "(→ liga); ↑C₁₂ engorda o caroço (→ empurra o vale pra fora e muda r₀)."),
        ("figuras/fig9_solucoes_radiais.png",
         "Fig. 9 — as soluções u(r) nos três regimes físicos",
         "EIXO X: r em fm. EIXO Y: u(r) normalizada para colar na reta em R. "
         "TRÊS PAINÉIS = três sistemas: n-n (a=−18,5 fm), unitariedade (a=∞), "
         "dêuteron (a=+5,4 fm). CURVAS COLORIDAS: os 4 potenciais; PONTILHADA: "
         "a reta g₀ = 1−r/a; TRACEJADAS: soluções analíticas (poço, e mPT no "
         "painel do meio). O QUE PROVA: fora de ~2 fm as 4 curvas viram A MESMA "
         "reta — impossível distinguir os potenciais. No painel do dêuteron a "
         "reta cruza zero em r=5,4 fm (o nó = assinatura do estado ligado); no "
         "unitário a reta é horizontal. O LJ (vermelho) nasce 'atrasado' (u=0 "
         "no caroço), dispara no vale e alcança as outras. SE MUDAR: ↑v em "
         "qualquer um → intercepto anda pra esquerda (a diminui)."),
        ("figuras/fig10_a_vs_intensidade.png",
         "Fig. 10 — a/r₀ contra a intensidade da interação",
         "EIXO X: o 'botão de força' — v para poço/mPT/gauss (painel a), C₆ "
         "para o LJ (painel b). EIXO Y: a/r₀ (adimensional). CURVAS: cheias = "
         "nosso código; tracejadas pretas/cinzas = fórmulas exatas Eq. 80 e "
         "117 (caem em cima = validação). POR QUE EXPLODE: cada divergência é "
         "um estado ligado nascendo — vindo da esquerda a→−∞, reaparece em "
         "+∞. Unitariedade = parar exatamente na parede (poço: v=1,2337; mPT: "
         "v=1; gauss: v=1,3420 — este último validado contra Jeszenszki et "
         "al. em 10⁻⁶). No LJ há VÁRIAS paredes: o poço de verdade comporta "
         "vários estados. SE MUDAR μ: as paredes NÃO andam (limiar só depende "
         "de v) — só a escala de a muda; é a invariância de escala."),
        ("figuras/fig_convergencia.png",
         "Convergência — o erro em função do passo Δr (log-log)",
         "EIXO X: passo da grade Δr em fm — PRA ESQUERDA = grade mais fina = "
         "conta mais cara. EIXO Y: erro relativo vs fórmula exata; 10⁻¹⁰ = 10 "
         "casas certas. LER ASSIM: reta no log-log = lei de potência; "
         "inclinação 2 → erro ∝ (Δr)², inclinação 4 → (Δr)⁴. PAINEL (a): no "
         "mPT suave o Numerov (quadrados) esmaga a diferença central — 10⁻¹⁰ "
         "contra 10⁻⁵. A SURPRESA: no poço descontínuo (azul) o Numerov vira "
         "reta de inclinação 1 (pior!) — a borda abrupta viola a hipótese de "
         "suavidade das correções dele. PAINEL (b): r₀ com trapézio vs "
         "Simpson. SE MUDAR: suavize a borda do poço e o Numerov volta a "
         "ganhar."),
        ("interativo/instantaneo_lab.png",
         "O laboratório interativo (laboratorio.bat)",
         "À esquerda: escolha do potencial e do método + sliders (intensidade, "
         "alcance, passo Δr em log10). NO MEIO: em cima V(r), embaixo u(r) com "
         "a reta 1−r/a. À DIREITA: a, r₀ (pelas duas quadraturas) e nº de nós, "
         "ao vivo. MELHOR EXPERIMENTO: arraste v devagar através do limiar e "
         "veja a reta deitar (a→∞) e o valor de a trocar de sinal — é a física "
         "da Fig. 10 nos seus dedos. Botões 'alvo:' rodam o ajuste automático."),
    ]
    corpo = passos + "".join(_img(*f) for f in figs)
    corpo += """
<h2>Os números que reproduzimos (e por que importam)</h2>
<table>
<tr><th>Resultado</th><th>Nosso valor</th><th>Referência</th><th>Leitura física</th></tr>
<tr><td>a (nêutron-nêutron)</td><td>−18,50 fm</td><td>−18,5 fm</td><td>quase liga: por pouco não existe o "dineutron"</td></tr>
<tr><td>a (dêuteron)</td><td>5,40 fm</td><td>5,4 fm</td><td>estado ligado raso, 1 nó em r=a</td></tr>
<tr><td>r₀ (unitariedade)</td><td>1,0000 fm</td><td>1,0 fm</td><td>a diverge mas r₀ continua finito e medível</td></tr>
<tr><td>E do dêuteron via (a, r₀)</td><td>−2,223 MeV</td><td>−2,224 MeV (exp.)</td><td>DOIS números preveem a energia com 0,05% de erro</td></tr>
<tr><td>E do dímero de ⁴He</td><td>−1,63 mK</td><td>−1,62 mK</td><td>mesma teoria, 9 ordens de grandeza de diferença em energia!</td></tr>
<tr><td>He₂ com LJ clássico</td><td>a = −178 Å (não liga)</td><td>ab initio: +90,4 Å (liga)</td><td>nosso achado: o dímero de hélio fica no fio da navalha</td></tr>
<tr><td>He₂ com Aziz HFD-B</td><td>a = +88,4 Å (liga!, 1 nó)</td><td>~88,5 Å (lit.)</td><td>o potencial realista resolve a navalha — E = −1,69 mK</td></tr>
<tr><td>n-p singleto (¹S₀)</td><td>estado virtual, −66 keV</td><td>a=−23,74 fm, r₀=2,77 fm</td><td>não existe "dêuteron singleto": κ&lt;0</td></tr>
<tr><td>³⁹K: 1º trímero de Efimov</td><td>B ≈ 403,5 G (previsto c/ a₁₋=−1500 a₀)</td><td>Zaccanti 2009</td><td>nosso módulo Feshbach + dado medido → o botão do laboratório</td></tr>
</table>
<p>Validação completa (4 camadas, 20 testes automáticos): aba <b>Resultados</b>.</p>"""
    return corpo


def gerar():
    abas = [
        ("comece", "🚀 Comece aqui", comece_aqui()),
        ("resultados", "📊 Resultados", _md("resultados/resumo.md")),
        ("guia", "📖 Guia (física)", _md("GUIA.md")),
        ("repo", "🗂️ Repositório", _md("GUIA_REPOSITORIO.md")),
        ("refs", "📚 Referências", _md("referencias/triagem_referencias.md")),
        ("visao", "🔭 Visão", _md("VISAO_PLATAFORMA.md")),
    ]
    botoes = "".join(
        f'<button class="aba" data-alvo="{i}">{titulo}</button>'
        for i, (id_, titulo, _) in enumerate(abas))
    paineis = "".join(
        f'<section class="painel" id="p{i}">{corpo}</section>'
        for i, (_, _, corpo) in enumerate(abas))

    html = f"""<!DOCTYPE html>
<html lang="pt-BR"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Laboratório de Espalhamento — Início</title>
<style>
:root{{--bg:#0e1117;--card:#161b22;--bd:#2d333f;--tx:#e6e8ee;--tx2:#9aa3b2;--ac:#58a6ff}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--tx);font:15px/1.6 -apple-system,"Segoe UI",Roboto,sans-serif}}
header{{padding:22px 24px 10px;max-width:960px;margin:0 auto}}
h1{{font-size:1.45rem;margin:0 0 4px}}
header p{{color:var(--tx2);margin:0 0 12px;font-size:.92rem}}
nav{{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--bd);padding:8px 24px;display:flex;gap:6px;flex-wrap:wrap;max-width:960px;margin:0 auto;z-index:5}}
.aba{{background:var(--card);color:var(--tx);border:1px solid var(--bd);border-radius:18px;padding:6px 13px;cursor:pointer;font-size:.86rem}}
.aba.on{{border-color:var(--ac);color:var(--ac)}}
main{{max-width:960px;margin:0 auto;padding:18px 24px 80px}}
.painel{{display:none}}.painel.on{{display:block}}
h2{{font-size:1.15rem;border-bottom:1px solid var(--bd);padding-bottom:6px;margin-top:28px}}
h3{{font-size:1rem}}
a{{color:var(--ac)}} code{{background:#0009;border-radius:4px;padding:1px 6px}}
pre{{background:#0d1420;border:1px solid var(--bd);border-radius:8px;padding:12px;overflow:auto}}
table{{border-collapse:collapse;width:100%;font-size:.88rem;margin:12px 0}}
th,td{{border:1px solid var(--bd);padding:7px 9px;text-align:left}}
th{{background:var(--card)}}
figure{{margin:26px 0;background:var(--card);border:1px solid var(--bd);border-radius:12px;padding:14px}}
figure img{{width:100%;border-radius:8px;background:#fff}}
figcaption{{margin-top:10px;font-size:.9rem;color:var(--tx2)}}
figcaption b{{color:var(--tx)}}
.passos li{{margin-bottom:10px}}
blockquote{{border-left:3px solid var(--ac);margin:0;padding:2px 14px;color:var(--tx2)}}
</style></head><body>
<header><h1>🧊 Laboratório de Espalhamento — tudo num lugar só</h1>
<p>Estudo do zero → resultados → guias → visão. Interativo: dois cliques em
<code>estudar.bat</code> (Jupyter) ou <code>laboratorio.bat</code> (sliders).
Linha do tempo: <a href="ROADMAP.html">ROADMAP.html</a>.</p></header>
<nav>{botoes}</nav>
<main>{paineis}</main>
<script>
const bs=[...document.querySelectorAll('.aba')],ps=[...document.querySelectorAll('.painel')];
function abre(i){{bs.forEach((b,j)=>b.classList.toggle('on',i==j));
ps.forEach((p,j)=>p.classList.toggle('on',i==j));window.scrollTo(0,0);
try{{localStorage.setItem('inicio_aba',i)}}catch(e){{}}}}
bs.forEach((b,i)=>b.onclick=()=>abre(i));
let ini=0;try{{ini=+(localStorage.getItem('inicio_aba')||0)}}catch(e){{}}
abre(ini);
</script></body></html>"""
    destino = os.path.join(AQUI, "INICIO.html")
    with open(destino, "w", encoding="utf-8") as f:
        f.write(html)
    print("  ->", destino, f"({os.path.getsize(destino)//1024} kB)")


if __name__ == "__main__":
    gerar()

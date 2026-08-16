# O Mac de simulação, do zero — checklist definitivo

> Objetivo: MacBook Intel i5 virando servidor de simulação estável, que
> você acessa de Sanca (ou de qualquer lugar) e onde o Claude trabalha
> direto no terminal. Instalação: ~40 min, a maior parte esperando.

## A máquina (identificada pela etiqueta)

**A1465 + EMC 2924 = MacBook Air 11", Early 2015.** O que isso muda:

| Peça | O que é | Situação no Linux |
|---|---|---|
| CPU | Core i5-5250U (Broadwell, 2 núcleos / 4 threads) | perfeito |
| Vídeo | Intel **HD 6000** (não é a 4000) | driver `i915` no kernel, sólido |
| RAM | 4 GB LPDDR3 **soldada** (não dá pra aumentar) | é o gargalo real |
| Wi-Fi | Broadcom BCM4360 (FCC QDS-BRCM1072) | **precisa do driver proprietário `wl`** |
| Firmware | EFI de 64 bits (2015) | ISO padrão boota liso |
| SSD | PCIe Apple (AHCI) | reconhecido normalmente |

Três consequências práticas, e elas mudam o roteiro:

1. **Boa notícia sobre a tela roxa:** Broadwell + HD 6000 é dos gráficos
   mais bem suportados que existem no Linux. O risco que te pegou no Mint
   cai muito — e o LTS derruba o resto.
2. **Este Mac NÃO tem porta de rede.** Só 2× USB 3.0 e Thunderbolt 2.
   Então "cabo de rede" não é uma opção sem adaptador — o plano A vira
   **celular no USB** (ver Parte 2).
3. **4 GB de RAM soldada** → marque **"Minimal installation"** no
   instalador. Sem isso o desktop come RAM que você quer pra simulação.

## Por que Ubuntu 24.04 LTS (e não Mint de novo)

O Mint é ótimo, mas pra máquina DEFINITIVA de simulação o Ubuntu
**24.04 LTS** ganha por três motivos práticos:

1. **LTS = atualizações conservadoras por 5 anos** — exatamente o
   antídoto do "apliquei update e deu tela roxa" (que foi quase
   certamente um kernel novo brigando com o vídeo do Mac).
2. É o padrão da computação científica: qualquer tutorial, pacote ou
   cluster (o Heaviside incluso) assume Ubuntu.
3. Melhor suporte testado a hardware de Mac (Wi-Fi Broadcom etc.).

## Parte 1 — preparar o pendrive (no seu desktop, 10 min)

1. Baixe a ISO: ubuntu.com/download/desktop → **Ubuntu 24.04 LTS**.
2. Abra o **balenaEtcher** (você já tem instalado!) → escolhe a ISO →
   escolhe o pendrive (8 GB+) → Flash.

## Parte 2 — instalar no Mac (30 min)

1. Espeta o pendrive no Mac → liga **segurando a tecla Option (Alt)** →
   escolhe **"EFI Boot"**.
2. ⚠️ Se a tela ficar bugada/roxa/preta no menu do instalador: escolha a
   opção **"Safe graphics"** (é pra isso que ela existe).
3. **Internet durante a instalação — o ovo e a galinha deste Mac.**
   O Wi-Fi BCM4360 só funciona depois que o driver proprietário entra, e
   o driver precisa de internet pra ser baixado. Como este Air não tem
   porta de rede:
   - **Plano A (o mais fácil):** Android no cabo USB → Configurações →
     *Ponto de acesso* → **"Ancoragem USB" / "USB tethering"**. O Ubuntu
     reconhece na hora, sem driver nenhum. iPhone também funciona.
   - **Plano B:** adaptador USB→Ethernet (qualquer um de R$ 40).
4. Marque **"Install third-party software for graphics and Wi-Fi
   hardware"** — é ELA que instala o `bcmwl-kernel-source`, o driver do
   seu Wi-Fi. Sem essa caixinha, o Mac fica mudo depois de instalado.
5. Marque **"Minimal installation"** — 4 GB de RAM não sobram pra
   LibreOffice e jogos que você nunca vai abrir.
6. "Erase disk and install Ubuntu" (apaga tudo) → segue o fluxo →
   usuário `pedro`, senha que você não esqueça.
7. Reiniciou e logou? Parte 3.

### Se o Wi-Fi não subir depois de instalado
Reconecta o celular no USB e roda:
```bash
sudo apt update && sudo apt install -y bcmwl-kernel-source
sudo modprobe wl
```
Se reclamar de Secure Boot, é só desabilitar no MOK (o Mac 2015 já vem
com Secure Boot desligado, então normalmente não aparece).

⚠️ Se DEPOIS de instalado a tela roxa voltar em algum boot: segure Shift
no boot → Advanced options → kernel anterior. E me chama, que aí ajusto
o nomodeset de forma permanente — é conserto de 2 minutos.

## Parte 3 — UMA linha e acabou a configuração

Abra o Terminal (Ctrl+Alt+T) e cole:

```bash
wget -qO- https://raw.githubusercontent.com/pemodest0/Low-energy-scattering/main/infra/setup_mac.sh | bash
```

*(enquanto o script não estiver no GitHub, copie o arquivo
`setup_mac.sh` desta pasta pro Mac num pendrive e rode `bash setup_mac.sh`)*

O script faz sozinho: atualiza o sistema, instala Python/git/SSH,
baixa o SEU laboratório do GitHub, roda os 25 testes pra provar que a
máquina calcula certo, e instala as duas peças da mágica:

- **Tailscale** → o Mac ganha um endereço fixo seguro, acessível de
  Sanca, do desktop, de onde for (rede privada própria, grátis).
- **Claude Code** → o Claude DENTRO do Mac: você abre o terminal
  (mesmo por SSH, de longe) digita `claude` e pede as coisas em
  português. Instalar pacote, rodar simulação, consertar driver — ele
  faz. É o "LLM cuidando do Linux" que você queria, só que sem gastar
  sua GPU: o cérebro fica na nuvem, as mãos ficam no Mac.

## Parte 4 — acessar de longe (5 min, uma vez)

1. No Mac: `sudo tailscale up` → abre um link → loga com Google.
2. No seu desktop: instala Tailscale também e loga na MESMA conta.
3. Pronto: de Sanca, `ssh pedro@mac-lab` e você está dentro. Digita
   `claude` e EU estou dentro. Simulação pesada? `heavyside` continua
   existindo pros trabalhos grandes; o Mac vira sua bancada 24h.

## Parte 5 — o detalhe que quebra o plano todo: a TAMPA

Notebook não foi feito pra ser servidor. Se você fechar a tampa, o
Ubuntu **suspende a máquina** — e o SSH de Sanca morre junto. O
`setup_mac.sh` já resolve isso (ignora a tampa e desliga suspensão),
mas é a linha mais importante do arquivo. Pra conferir depois:

```bash
grep HandleLidSwitch /etc/systemd/logind.conf   # deve dizer =ignore
systemctl status sleep.target                   # deve dizer "masked"
```

Mantenha na tomada. Bateria de 2015 em ciclo eterno de 100% envelhece
rápido, mas é peça de consumo — o que você não quer é a máquina
desligando sozinha quando você está a 200 km.

## Expectativa honesta de desempenho

O i5-5250U tem 2 núcleos. Ele **não** é o músculo — o Heaviside é. O
que este Mac faz muito bem: rodar os 25 testes, varreduras leves do
laboratório, o app Streamlit, e ficar de pé 24h esperando você pedir
coisas pelo `claude`. Se um cálculo passar de ~20 min aqui, é sinal de
que ele pertence ao cluster. E com 4 GB, cuidado com grades grandes:
o script liga **zram** (compressão de RAM) justamente pra dar fôlego.

## Sobre a LLM local no desktop

Adio com carinho: um i5 antigo não roda nada útil, e no desktop uma
LLM local pequena faria pior o que o Claude Code já faz melhor (e o
seu negócio é a física, não babá de modelo). Se um dia quiser
brincar com isso por hobby, a gente monta — mas pro laboratório, o
caminho acima é o profissional.

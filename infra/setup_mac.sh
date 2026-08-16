#!/usr/bin/env bash
# ============================================================
#  setup_mac.sh — transforma um Ubuntu recém-instalado no
#  servidor de simulação do laboratório. Rode UMA vez:
#      bash setup_mac.sh
#  Idempotente: rodar de novo não estraga nada.
#  Alvo: MacBook Air 11" Early 2015 (A1465 / EMC 2924)
#        i5-5250U Broadwell, HD 6000, 4 GB RAM, Wi-Fi BCM4360
# ============================================================
set -e
echo "==> [1/8] Atualizando o sistema..."
sudo apt update && sudo apt upgrade -y

echo "==> [2/8] Ferramentas essenciais (python, git, ssh, utilitarios)..."
sudo apt install -y git python3-pip python3-venv build-essential \
    openssh-server curl wget htop tmux
sudo systemctl enable --now ssh

echo "==> [3/8] Driver Wi-Fi Broadcom BCM4360 (o Wi-Fi deste Mac)..."
sudo apt install -y bcmwl-kernel-source || echo "   (sem Broadcom, ok)"

echo "==> [4/8] Modo servidor: ignorar tampa fechada e nunca suspender..."
# Sem isto, fechar a tampa derruba o SSH que voce usa de Sanca.
sudo sed -i 's/^#*HandleLidSwitch=.*/HandleLidSwitch=ignore/'                 /etc/systemd/logind.conf
sudo sed -i 's/^#*HandleLidSwitchExternalPower=.*/HandleLidSwitchExternalPower=ignore/' /etc/systemd/logind.conf
sudo sed -i 's/^#*HandleLidSwitchDocked=.*/HandleLidSwitchDocked=ignore/'     /etc/systemd/logind.conf
grep -q '^HandleLidSwitch=' /etc/systemd/logind.conf || \
    echo -e "HandleLidSwitch=ignore\nHandleLidSwitchExternalPower=ignore\nHandleLidSwitchDocked=ignore" \
    | sudo tee -a /etc/systemd/logind.conf >/dev/null
sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target
sudo systemctl restart systemd-logind || true

echo "==> [5/8] zram: folego de RAM para os 4 GB soldados deste Air..."
sudo apt install -y zram-tools
echo -e "ALGO=zstd\nPERCENT=60" | sudo tee /etc/default/zramswap >/dev/null
sudo systemctl enable --now zramswap.service || true

echo "==> [6/8] Tailscale (acesso remoto seguro de qualquer lugar)..."
curl -fsSL https://tailscale.com/install.sh | sh

echo "==> [7/8] Node + Claude Code (o Claude dentro desta maquina)..."
curl -fsSL https://deb.nodesource.com/setup_22.x | sudo -E bash -
sudo apt install -y nodejs
sudo npm install -g @anthropic-ai/claude-code

echo "==> [8/8] O laboratorio (do SEU GitHub) + prova de fogo..."
if [ ! -d "$HOME/lab" ]; then
    git clone https://github.com/pemodest0/Low-energy-scattering.git "$HOME/lab"
fi
cd "$HOME/lab"
python3 -m venv .venv
source .venv/bin/activate
pip install --quiet -r requirements.txt pytest
python -m pytest tests/ -q

echo ""
echo "============================================================"
echo "  MAQUINA PRONTA. Os 25 testes acima devem estar verdes."
echo "  Proximos passos (uma vez só):"
echo "    1) sudo tailscale up     (login -> acesso remoto)"
echo "    2) claude                (login -> o Claude no terminal)"
echo "  Depois disso, de Sanca:  ssh $(whoami)@<nome-no-tailscale>"
echo "============================================================"

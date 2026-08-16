#!/usr/bin/env bash
# =============================================================================
#  setup_lab.sh — instala o Laboratório de Espalhamento a Baixas Energias
#
#  Alvos suportados:
#    * MacBook Intel rodando Ubuntu 24.04 LTS   -> bancada 24h ("mac")
#    * WSL2 Ubuntu no desktop Windows           -> espelho local ("wsl")
#    * Qualquer Ubuntu/Debian                   -> genérico ("linux")
#
#  Uso:
#    bash setup_lab.sh                 # detecta o alvo sozinho
#    bash setup_lab.sh --role=mac      # força o alvo
#    bash setup_lab.sh --no-tailscale  # pula a rede privada
#    bash setup_lab.sh --repo-only     # só clona e testa o laboratório
#
#  É idempotente: pode rodar de novo quantas vezes quiser.
# =============================================================================

set -Eeuo pipefail

REPO_URL="https://github.com/pemodest0/Low-energy-scattering.git"
LAB_DIR="${LAB_DIR:-$HOME/lab/Low-energy-scattering}"
VENV_DIR="$LAB_DIR/.venv"
LOG="$HOME/setup_lab.log"

ROLE=""
DO_TAILSCALE=1
REPO_ONLY=0

# ---------------------------------------------------------------- aparência --
c_ok()   { printf '\033[1;32m  ok  \033[0m %s\n' "$*"; }
c_info() { printf '\033[1;36m ---- \033[0m %s\n' "$*"; }
c_warn() { printf '\033[1;33m aten \033[0m %s\n' "$*"; }
c_err()  { printf '\033[1;31m ERRO \033[0m %s\n' "$*" >&2; }
step()   { printf '\n\033[1;35m==> %s\033[0m\n' "$*"; }

trap 'c_err "falhou na linha $LINENO. Log completo em $LOG"' ERR

# ---------------------------------------------------------------- argumentos --
for arg in "$@"; do
  case "$arg" in
    --role=*)       ROLE="${arg#*=}" ;;
    --no-tailscale) DO_TAILSCALE=0 ;;
    --repo-only)    REPO_ONLY=1 ;;
    -h|--help)      sed -n '2,20p' "$0"; exit 0 ;;
    *)              c_warn "argumento ignorado: $arg" ;;
  esac
done

exec > >(tee -a "$LOG") 2>&1
echo "===== setup_lab.sh — $(date -Is) ====="

# ---------------------------------------------------------------- pré-checks --
[[ $EUID -eq 0 ]] && { c_err "não rode como root. Rode como você mesmo; o script pede sudo quando precisa."; exit 1; }
command -v apt-get >/dev/null || { c_err "isto não é um Ubuntu/Debian."; exit 1; }

# Mantém o sudo acordado durante a instalação inteira.
sudo -v
( while true; do sudo -n true; sleep 50; kill -0 "$$" 2>/dev/null || exit; done ) &
SUDO_KEEPALIVE=$!
trap 'kill $SUDO_KEEPALIVE 2>/dev/null || true' EXIT

# ---------------------------------------------------------------- detecção ---
detect_role() {
  if grep -qi microsoft /proc/version 2>/dev/null; then echo wsl; return; fi
  if [[ -d /sys/devices/platform/applesmc.768 ]] \
     || sudo dmidecode -s system-manufacturer 2>/dev/null | grep -qi apple; then
    echo mac; return
  fi
  echo linux
}
[[ -z "$ROLE" ]] && ROLE="$(detect_role)"
c_info "alvo: $ROLE   |   laboratório em: $LAB_DIR"

# =============================================================================
step "1/9  Sistema base"
# =============================================================================
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends \
  build-essential gfortran pkg-config \
  git curl wget ca-certificates gnupg \
  python3 python3-venv python3-dev python3-pip \
  libopenblas-dev liblapack-dev \
  tmux htop rsync unzip jq \
  fonts-dejavu-core
c_ok "pacotes base instalados"

# =============================================================================
step "2/9  Ajustes específicos do Mac"
# =============================================================================
if [[ "$ROLE" == "mac" ]]; then

  # --- Wi-Fi Broadcom: o passo que 90% das pessoas esquece -------------------
  if lspci -nn 2>/dev/null | grep -qi 'network.*broadcom'; then
    if ! lsmod | grep -q '^wl'; then
      c_info "Wi-Fi Broadcom detectado sem driver — instalando bcmwl"
      sudo apt-get install -y bcmwl-kernel-source || \
        c_warn "bcmwl falhou; use cabo/USB e veja o CHECKLIST_MAC.md"
    else
      c_ok "driver Wi-Fi Broadcom já carregado"
    fi
  fi

  # --- Ventoinha: Mac com Linux não controla a ventoinha sozinho ------------
  # Só o mbpfan. Instalar mbpfan e macfanctld juntos faz os dois brigarem
  # pelo mesmo controle e a ventoinha oscilar.
  sudo apt-get install -y mbpfan lm-sensors || c_warn "mbpfan indisponível"
  if systemctl list-unit-files | grep -q '^mbpfan'; then
    sudo systemctl enable --now mbpfan || true
    c_ok "controle de ventoinha ativo (a máquina não vai cozinhar)"
  fi

  # --- Servidor de tampa fechada: NÃO dormir --------------------------------
  sudo mkdir -p /etc/systemd/logind.conf.d
  sudo tee /etc/systemd/logind.conf.d/99-lab-servidor.conf >/dev/null <<'EOF'
# A bancada fica de tampa fechada e não dorme. Ela é servidor.
[Login]
HandleLidSwitch=ignore
HandleLidSwitchDocked=ignore
HandleLidSwitchExternalPower=ignore
IdleAction=ignore
EOF
  sudo systemctl restart systemd-logind || true
  sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target 2>/dev/null || true
  c_ok "suspensão desligada (tampa fechada = continua rodando)"

  # --- Rede de segurança contra a tela roxa ---------------------------------
  # O menu do GRUB visível é o seguro: se um kernel novo brigar com o vídeo,
  # você reinicia, escolhe "Advanced options" e volta pro kernel anterior.
  if [[ -f /etc/default/grub ]]; then
    sudo sed -i 's/^GRUB_TIMEOUT_STYLE=.*/GRUB_TIMEOUT_STYLE=menu/' /etc/default/grub
    grep -q '^GRUB_TIMEOUT_STYLE=' /etc/default/grub || \
      echo 'GRUB_TIMEOUT_STYLE=menu' | sudo tee -a /etc/default/grub >/dev/null
    sudo sed -i 's/^GRUB_TIMEOUT=.*/GRUB_TIMEOUT=5/' /etc/default/grub
    sudo update-grub || c_warn "update-grub falhou (não é fatal)"
    c_ok "menu do GRUB visível por 5 s — sua saída de emergência da tela roxa"
  fi

  # Guarda os últimos 3 kernels em vez de podar tudo.
  echo 'APT::NeverAutoRemove:: "^linux-image-.*";' | \
    sudo tee /etc/apt/apt.conf.d/99-keep-kernels >/dev/null
  c_ok "kernels antigos preservados para rollback"

elif [[ "$ROLE" == "wsl" ]]; then
  c_info "WSL2: pulando ventoinha, tampa e GRUB (o Windows cuida disso)"
  # Torna o laboratório visível do Windows em \\wsl$\Ubuntu\home\...
  c_ok "nada a ajustar"
else
  c_info "alvo genérico: pulando ajustes de Mac"
fi

# =============================================================================
step "3/9  Laboratório: clone + ambiente Python"
# =============================================================================
mkdir -p "$(dirname "$LAB_DIR")"
if [[ -d "$LAB_DIR/.git" ]]; then
  c_info "repositório já existe — atualizando"
  git -C "$LAB_DIR" pull --ff-only || c_warn "pull não fez fast-forward; resolva à mão depois"
else
  git clone --depth 1 "$REPO_URL" "$LAB_DIR"
fi
c_ok "repositório em $LAB_DIR"

python3 -m venv "$VENV_DIR" 2>/dev/null || true
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip wheel setuptools -q
if [[ -f "$LAB_DIR/requirements.txt" ]]; then
  python -m pip install -r "$LAB_DIR/requirements.txt" -q
else
  c_warn "requirements.txt não encontrado — instalando o mínimo"
  python -m pip install numpy scipy matplotlib pandas pytest pyyaml streamlit jupyterlab -q
fi
python -m pip install pytest jupyterlab -q
c_ok "ambiente Python pronto ($(python --version))"

# =============================================================================
step "4/9  Prova de fogo: os 25 testes"
# =============================================================================
cd "$LAB_DIR"
if python -m pytest tests/ -q; then
  c_ok "todos os testes passaram — a máquina reproduz a física"
else
  c_warn "algum teste falhou. A instalação continua, mas ANOTE isso."
fi

[[ $REPO_ONLY -eq 1 ]] && { c_ok "modo --repo-only: terminado."; exit 0; }

# =============================================================================
step "5/9  Acesso remoto: SSH"
# =============================================================================
sudo apt-get install -y openssh-server
sudo systemctl enable --now ssh
# Só chave, sem senha, é o certo — mas só depois que você tiver copiado a sua.
if [[ -f "$HOME/.ssh/authorized_keys" ]]; then
  sudo sed -i 's/^#\?PasswordAuthentication.*/PasswordAuthentication no/' /etc/ssh/sshd_config
  sudo systemctl reload ssh
  c_ok "SSH por chave apenas (senha desligada)"
else
  c_warn "SSH ativo COM senha — copie sua chave (ssh-copy-id) e rode o script de novo"
fi

# =============================================================================
step "6/9  Rede privada: Tailscale"
# =============================================================================
if [[ $DO_TAILSCALE -eq 1 ]]; then
  if ! command -v tailscale >/dev/null; then
    curl -fsSL https://tailscale.com/install.sh | sh
  fi
  sudo systemctl enable --now tailscaled || true
  c_ok "Tailscale instalado"
  c_warn "FALTA VOCÊ:  sudo tailscale up   (abre um link, você loga uma vez)"
else
  c_info "Tailscale pulado por opção"
fi

# =============================================================================
step "7/9  Claude Code no terminal"
# =============================================================================
if ! command -v node >/dev/null || [[ "$(node -v 2>/dev/null | tr -d 'v' | cut -d. -f1)" -lt 18 ]]; then
  curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
  sudo apt-get install -y nodejs
fi
mkdir -p "$HOME/.npm-global"
npm config set prefix "$HOME/.npm-global"
grep -q 'npm-global/bin' "$HOME/.bashrc" || \
  echo 'export PATH="$HOME/.npm-global/bin:$PATH"' >> "$HOME/.bashrc"
export PATH="$HOME/.npm-global/bin:$PATH"
npm install -g @anthropic-ai/claude-code || c_warn "instalação do Claude Code falhou — tente à mão depois"
command -v claude >/dev/null && c_ok "Claude Code instalado ($(claude --version 2>/dev/null || echo ok))"

# =============================================================================
step "8/9  Serviços do laboratório (JupyterLab + Streamlit)"
# =============================================================================
if [[ "$ROLE" != "wsl" ]]; then
  sudo tee /etc/systemd/system/lab-jupyter.service >/dev/null <<EOF
[Unit]
Description=JupyterLab do laboratorio de espalhamento
After=network-online.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$LAB_DIR
Environment="PATH=$VENV_DIR/bin:/usr/local/bin:/usr/bin:/bin"
ExecStart=$VENV_DIR/bin/jupyter lab --no-browser --ip=127.0.0.1 --port=8888
Restart=on-failure
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

  if [[ -f "$LAB_DIR/app/Inicio.py" ]]; then
    sudo tee /etc/systemd/system/lab-streamlit.service >/dev/null <<EOF
[Unit]
Description=App Streamlit do laboratorio de espalhamento
After=network-online.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$LAB_DIR
Environment="PATH=$VENV_DIR/bin:/usr/local/bin:/usr/bin:/bin"
ExecStart=$VENV_DIR/bin/streamlit run app/Inicio.py --server.address=0.0.0.0 --server.port=8501 --server.headless=true
Restart=on-failure
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF
  fi

  sudo systemctl daemon-reload
  sudo systemctl enable --now lab-jupyter || c_warn "lab-jupyter não subiu"
  [[ -f /etc/systemd/system/lab-streamlit.service ]] && \
    { sudo systemctl enable --now lab-streamlit || c_warn "lab-streamlit não subiu"; }
  c_ok "serviços criados (sobem sozinhos a cada boot)"
else
  c_info "WSL2 não tem systemd por padrão — use os atalhos 'lab-app' e 'lab-jup'"
fi

# =============================================================================
step "9/9  Atalhos de conveniência"
# =============================================================================
MARK="# >>> laboratorio de espalhamento >>>"
if ! grep -q "$MARK" "$HOME/.bashrc"; then
cat >> "$HOME/.bashrc" <<EOF

$MARK
export LAB_DIR="$LAB_DIR"
alias lab='cd "\$LAB_DIR" && source .venv/bin/activate'
alias lab-test='(cd "\$LAB_DIR" && .venv/bin/python -m pytest tests/ -q)'
alias lab-run='(cd "\$LAB_DIR" && .venv/bin/python main.py --sem-ajustes)'
alias lab-app='(cd "\$LAB_DIR" && .venv/bin/streamlit run app/Inicio.py)'
alias lab-jup='(cd "\$LAB_DIR" && .venv/bin/jupyter lab)'
alias lab-up='(cd "\$LAB_DIR" && git pull --ff-only && .venv/bin/pip install -r requirements.txt -q && .venv/bin/python -m pytest tests/ -q)'
# <<< laboratorio de espalhamento <<<
EOF
fi
c_ok "atalhos: lab, lab-test, lab-run, lab-app, lab-jup, lab-up"

# =============================================================================
echo
printf '\033[1;32m'
cat <<'EOF'
=============================================================================
  Laboratório instalado.
=============================================================================
EOF
printf '\033[0m'

cat <<EOF

  Laboratório .... $LAB_DIR
  Log ............ $LOG
  Recarregue o shell:  source ~/.bashrc

  Falta você fazer, uma vez só:
    1) sudo tailscale up          -> abre um link, loga, e a máquina entra na sua rede
    2) claude                     -> loga na sua conta Anthropic
    3) tailscale ip -4            -> anote o IP; é por ele que você entra de Sanca

  Depois disso, de qualquer lugar:
    ssh $USER@<ip-tailscale>
    claude

EOF

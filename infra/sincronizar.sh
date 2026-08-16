#!/usr/bin/env bash
# =============================================================================
#  sincronizar.sh  —  Mac / Ubuntu
#
#  Um comando so, nos dois sentidos:
#     1. baixa o que a outra maquina fez
#     2. sobe o que esta maquina fez
#
#  Uso:   ./infra/sincronizar.sh  "mensagem opcional"
#  Regra: rode ANTES de comecar a trabalhar e DEPOIS de terminar.
# =============================================================================
set -u
cd "$(dirname "$0")/.." || exit 1
MSG="${1:-sincronizacao automatica de $(hostname -s) em $(date '+%d/%m/%Y %H:%M')}"

echo "============================================================"
echo "  SINCRONIZAR   pasta: $(pwd)"
echo "  maquina: $(hostname -s)"
echo "============================================================"

command -v git >/dev/null || { echo "[ERRO] git nao instalado."; exit 1; }
[ -n "$(git config user.name || true)" ] || {
  git config --global user.name  "Pedro Henrique G. Modesto"
  git config --global user.email "pedrohenriquemodesto4@gmail.com"; }

# ---- 1. guarda o que esta maquina fez, ANTES de baixar -----------------------
git add -A
if ! git diff --cached --quiet; then
    git commit -m "$MSG" -q && echo "[ok] mudancas locais guardadas no commit"
else
    echo "[--] nada novo nesta maquina"
fi

# ---- 2. baixa o que a outra maquina fez -------------------------------------
echo "[--] baixando do GitHub..."
if ! git pull --rebase origin main; then
    echo
    echo "############################################################"
    echo "  CONFLITO: as duas maquinas mexeram no MESMO arquivo."
    echo
    echo "  Veja quais:   git status"
    echo "  Abra o arquivo, escolha o que fica, apague as marcas"
    echo "  <<<<<<<  =======  >>>>>>>  e entao:"
    echo "     git add <arquivo>"
    echo "     git rebase --continue"
    echo
    echo "  Para desistir e voltar ao estado de antes:"
    echo "     git rebase --abort"
    echo "############################################################"
    exit 1
fi
echo "[ok] atualizado com a outra maquina"

# ---- 3. sobe -----------------------------------------------------------------
git push origin main && echo "[ok] enviado ao GitHub" || { echo "[ERRO] push falhou"; exit 1; }

echo
echo "PRONTO. As duas maquinas estao iguais."
echo "  ultimo commit: $(git log -1 --pretty='%h  %s')"

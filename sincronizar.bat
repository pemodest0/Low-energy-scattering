@echo off
chcp 65001 >nul
setlocal EnableDelayedExpansion

REM ============================================================
REM  sincronizar.bat  -  Windows (desktop)
REM
REM  Um comando so, nos DOIS sentidos:
REM     1. guarda o que voce fez aqui
REM     2. baixa o que voce fez no Mac
REM     3. sobe tudo
REM
REM  DOIS CLIQUES neste arquivo.
REM  Regra: rode ANTES de comecar e DEPOIS de terminar.
REM
REM  No Mac o equivalente e:   ./infra/sincronizar.sh
REM ============================================================

cd /d "%~dp0"
set REPO=https://github.com/pemodest0/Low-energy-scattering.git

echo.
echo ============================================================
echo   SINCRONIZAR
echo   pasta:   %CD%
echo   maquina: %COMPUTERNAME%  (Windows)
echo ============================================================
echo.

REM ---------- git existe? ----------
where git >nul 2>&1
if errorlevel 1 (
    echo [ERRO] git nao encontrado. Instale: https://git-scm.com/download/win
    pause & exit /b 1
)

REM ---------- identidade ----------
for /f "delims=" %%N in ('git config --global user.name 2^>nul') do set GITNAME=%%N
if "!GITNAME!"=="" (
    git config --global user.name "Pedro Henrique G. Modesto"
    git config --global user.email "pedrohenriquemodesto4@gmail.com"
)

REM ---------- repositorio saudavel? ----------
git rev-parse --git-dir >nul 2>&1
if errorlevel 1 (
    echo [ERRO] O repositorio esta quebrado.
    echo   Rode publicar.bat primeiro: ele conserta sozinho.
    pause & exit /b 1
)

REM ---------- 1. guarda o que ESTA maquina fez ----------
echo ------------------------------------------------------------
echo   1. O QUE MUDOU AQUI
echo ------------------------------------------------------------
git add -A
git status --short
echo.

git diff --cached --quiet
if errorlevel 1 (
    set MSG=%~1
    if "!MSG!"=="" set MSG=sincronizacao do desktop %COMPUTERNAME% em %DATE% %TIME:~0,5%
    git commit -m "!MSG!" -q
    echo [ok] mudancas locais guardadas
) else (
    echo [--] nada novo nesta maquina
)

REM ---------- 2. baixa o que a OUTRA maquina fez ----------
echo.
echo ------------------------------------------------------------
echo   2. BAIXANDO O QUE VOCE FEZ NO MAC
echo ------------------------------------------------------------
git pull --rebase origin main
if errorlevel 1 (
    echo.
    echo ############################################################
    echo   CONFLITO: as duas maquinas mexeram no MESMO arquivo.
    echo.
    echo   Veja quais:   git status
    echo   Abra o arquivo, escolha o que fica, apague as marcas
    echo   ^<^<^<^<^<^<^<   =======   ^>^>^>^>^>^>^>   e entao:
    echo        git add ARQUIVO
    echo        git rebase --continue
    echo.
    echo   Para desistir e voltar ao estado anterior:
    echo        git rebase --abort
    echo ############################################################
    pause & exit /b 1
)
echo [ok] atualizado com o Mac

REM ---------- 3. sobe ----------
echo.
echo ------------------------------------------------------------
echo   3. ENVIANDO
echo ------------------------------------------------------------
git push origin main
if errorlevel 1 (
    echo [ERRO] push falhou. Mensagem acima.
    pause & exit /b 1
)

echo.
echo ============================================================
echo   PRONTO. As duas maquinas estao iguais.
echo.
git log -1 --pretty="   ultimo commit: %%h  %%s"
echo.
echo   No Mac, para pegar isto:
echo      cd ~/lab/Low-energy-scattering ^&^& ./infra/sincronizar.sh
echo ============================================================
echo.
pause

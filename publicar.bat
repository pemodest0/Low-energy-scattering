@echo off
chcp 65001 >nul
setlocal EnableDelayedExpansion

REM ============================================================
REM  publicar.bat  -  CONSERTA o git e sobe o laboratorio.
REM
REM  >>> No dia a dia use sincronizar.bat, nao este. <<<
REM  Este aqui e para quando o repositorio quebra: ele reconstroi
REM  a pasta .git do zero sem tocar nos seus arquivos.
REM
REM  DOIS CLIQUES neste arquivo. Ele acha a propria pasta.
REM
REM  Se a pasta .git estiver quebrada (foi o caso em 15/08/2026:
REM  faltava .git\objects), ele CONSERTA sozinho, guardando a
REM  antiga como .git_quebrado em vez de apagar.
REM ============================================================

set REPO=https://github.com/pemodest0/Low-energy-scattering.git

cd /d "%~dp0"

echo.
echo ============================================================
echo   PUBLICAR O LABORATORIO
echo   pasta: %CD%
echo ============================================================
echo.

REM ---------- 1. o git existe? ----------
where git >nul 2>&1
if errorlevel 1 (
    echo [ERRO] O comando "git" nao foi encontrado.
    echo   Instale em: https://git-scm.com/download/win
    pause
    exit /b 1
)
echo [ok] git encontrado

REM ---------- 2. quem e voce? ----------
for /f "delims=" %%N in ('git config --global user.name 2^>nul') do set GITNAME=%%N
if "!GITNAME!"=="" (
    echo [--] configurando seu nome e email pela primeira vez
    git config --global user.name "Pedro Henrique G. Modesto"
    git config --global user.email "pedrohenriquemodesto4@gmail.com"
)
echo [ok] identidade configurada

REM ---------- 3. o repositorio funciona? ----------
git rev-parse --git-dir >nul 2>&1
if errorlevel 1 goto CONSERTAR
echo [ok] repositorio saudavel
goto TEM_REMOTO

:CONSERTAR
echo.
echo ------------------------------------------------------------
echo   [!] O repositorio esta quebrado. Consertando.
echo ------------------------------------------------------------
if exist ".git" (
    if exist ".git_quebrado" rmdir /s /q ".git_quebrado"
    move ".git" ".git_quebrado" >nul 2>&1
    if exist ".git" (
        echo [ERRO] Nao consegui mover a pasta .git
        echo   Feche VS Code, GitHub Desktop e terminais, e rode de novo.
        pause
        exit /b 1
    )
    echo   [ok] .git antiga guardada como .git_quebrado
)

git init -b main >nul 2>&1
if errorlevel 1 ( git init >nul 2>&1 & git checkout -b main >nul 2>&1 )
echo   [ok] repositorio novo criado

git remote add origin %REPO% >nul 2>&1
echo   [ok] endereco do GitHub configurado

echo   [--] baixando o historico do GitHub (pode pedir login)...
git fetch origin main
if errorlevel 1 (
    echo.
    echo [ERRO] Nao consegui baixar do GitHub. Causas comuns:
    echo   - sem internet
    echo   - falta autenticar: na primeira vez o Git abre uma janela
    echo     do navegador. Se nao abriu, instale o Git Credential Manager
    echo     (vem junto com o Git for Windows).
    echo.
    pause
    exit /b 1
)
echo   [ok] historico baixado

git reset origin/main >nul 2>&1
echo   [ok] repositorio religado ao GitHub, sem tocar nos seus arquivos
echo.

:TEM_REMOTO
REM ---------- 4. tem endereco do GitHub? ----------
git remote get-url origin >nul 2>&1
if errorlevel 1 (
    git remote add origin %REPO%
    echo [ok] endereco do GitHub adicionado
)

REM ---------- 5. o que mudou ----------
echo.
echo ------------------------------------------------------------
echo   O QUE VAI SUBIR
echo ------------------------------------------------------------
git add notebooks src tests referencias notas_teoria figuras infra app resultados 2>nul
git add .github 2>nul
for %%F in (MAPA.md COMO_TRABALHAR_NAS_DUAS_MAQUINAS.md HISTORICO.md VISITA_PATRICIA.md Makefile LICENSE CITATION.cff requirements.txt .gitignore .gitattributes publicar.bat sincronizar.bat README.md main.py pyproject.toml) do (
    if exist "%%F" git add "%%F" 2>nul
)
git status --short
echo.

REM ---------- 6. tem o que commitar? ----------
git diff --cached --quiet
if not errorlevel 1 (
    echo [aviso] Nada novo para commitar.
    pause
    exit /b 0
)

REM ---------- 7. commit ----------
git commit -m "Laboratorio de 2 corpos guiado, trimero no hiper-raio, incerteza e portoes de validade"
if errorlevel 1 (
    echo [ERRO] O commit falhou. Mensagem acima.
    pause
    exit /b 1
)
echo [ok] commit feito

REM ---------- 8. push ----------
echo.
echo Enviando para o GitHub...
git push -u origin main
if errorlevel 1 (
    echo.
    echo [ERRO] O push falhou. Se disser que o remoto esta na frente:
    echo     git pull --rebase origin main
    echo   e rode este arquivo de novo.
    pause
    exit /b 1
)

echo.
echo ============================================================
echo   PRONTO.
echo   https://github.com/pemodest0/Low-energy-scattering
echo.
echo   No Mac:  cd ~/lab/Low-energy-scattering ^&^& git pull
echo.
echo   A pasta .git_quebrado pode ser apagada quando voce quiser.
echo ============================================================
echo.
pause

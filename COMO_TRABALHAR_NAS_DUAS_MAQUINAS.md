# Trabalhar no Mac e no desktop ao mesmo tempo

O objetivo: rodar uma simulação no Mac e ter o resultado no desktop, ou o
contrário, sem copiar nada em pendrive e sem nunca ficar em dúvida sobre qual
versão é a boa.

**A ponte é o GitHub.** Ele já existe, é gratuito, guarda o histórico inteiro e
funciona nos dois sistemas. Não precisa de nada novo.

---

## O ciclo, e é só isto

```
   MAC                    GitHub                  DESKTOP
    │                        │                        │
    │  ./infra/              │                        │
    │  sincronizar.sh   ───► │ ──►  sincronizar.bat   │
    │                        │      (dois cliques)    │
```

| Máquina | O que você faz |
|---|---|
| **Windows** | dois cliques em `sincronizar.bat` |
| **Mac / Ubuntu** | `cd ~/lab/Low-energy-scattering && ./infra/sincronizar.sh` |

Os dois fazem exatamente a mesma coisa, nesta ordem:

1. **guarda** o que você fez nesta máquina (commit)
2. **baixa** o que você fez na outra (pull com rebase)
3. **sobe** tudo (push)

---

## A única regra que importa

> **Sincronize ANTES de começar a trabalhar e DEPOIS de terminar.**

Se você fizer isso, nunca haverá conflito. Conflito só acontece quando as duas
máquinas mexem no **mesmo arquivo** sem sincronizar no meio.

Cole isso na cabeça como se fosse abrir e fechar uma porta: você abre (sincroniza)
ao entrar, trabalha, e fecha (sincroniza) ao sair.

---

## Primeira vez no Mac

```bash
mkdir -p ~/lab && cd ~/lab
git clone https://github.com/pemodest0/Low-energy-scattering.git
cd Low-energy-scattering
chmod +x infra/sincronizar.sh infra/setup_mac.sh
./infra/setup_mac.sh          # instala Python, numpy, scipy, matplotlib, pytest
python3 -m pytest tests/ -q   # confere que tudo passa
```

A partir daí, só `./infra/sincronizar.sh`.

---

## Rodar uma simulação longa no Mac e pegar aqui

O Mac é a bancada 24/7 — é para isso que ele existe. O fluxo:

**No Mac**, deixe rodando sem depender da janela do terminal ficar aberta:

```bash
cd ~/lab/Low-energy-scattering
nohup python3 -m src.relatorios.gerar_resultados > resultados/log_$(date +%F).txt 2>&1 &
```

Acompanhar depois:

```bash
tail -f resultados/log_*.txt
```

Quando terminar, ainda no Mac:

```bash
./infra/sincronizar.sh "resultados da varredura de alcance efetivo"
```

**No desktop**, dois cliques em `sincronizar.bat`. Os CSVs novos aparecem em
`resultados/`.

---

## O que sobe e o que não sobe

| Sobe | Não sobe |
|---|---|
| `src/`, `tests/`, `notebooks/` | `arquivo/` (coisas aposentadas, 12 MB) |
| `referencias/`, `notas_teoria/` | `__pycache__/`, `.pytest_cache/` |
| `figuras/`, `resultados/` (CSVs) | lixo de compilação do LaTeX (`.aux`, `.log`, `.fls`) |
| `MAPA.md`, `HISTORICO.md` | `.git_quebrado/` |

Está tudo no `.gitignore`.

**Se um dia você gerar um arquivo de resultado com centenas de megabytes**, não
deixe entrar no git — ele guarda **todas** as versões para sempre e o repositório
incha sem volta. Nesse caso me avise que a gente resolve (Git LFS ou pasta
separada). CSVs de alguns megabytes estão perfeitamente bem.

---

## Se der conflito

Acontece quando as duas máquinas mudaram a **mesma linha** do **mesmo arquivo**.
O script para e mostra o que fazer. Em resumo:

```bash
git status                    # quais arquivos conflitaram
```

Abra o arquivo. Você vai ver:

```
<<<<<<< HEAD
    o que ESTA maquina escreveu
=======
    o que a OUTRA maquina escreveu
>>>>>>> abc1234
```

Escolha o que fica, apague as três linhas de marcação, e:

```bash
git add ARQUIVO
git rebase --continue
```

Para desistir e voltar exatamente ao estado de antes, sem perder nada:

```bash
git rebase --abort
```

---

## Por que o `.gitattributes` existe

Windows termina linha com CRLF, Mac com LF. Sem normalizar, o git jura que você
mudou o arquivo inteiro quando você não mudou nada — e aí *tudo* vira conflito.
O `.gitattributes` na raiz resolve isso. Não mexa nele.

---

## Sobre os notebooks

Notebook (`.ipynb`) guarda a saída junto com o código. Se você rodar o mesmo
notebook nas duas máquinas, o git vê arquivos diferentes mesmo com o código igual.

**Recomendação simples:** antes de sincronizar, no Jupyter/VS Code faça
*Kernel → Restart & Clear All Outputs* no notebook que você mexeu. Aí o diff fica
só do código, que é o que interessa.

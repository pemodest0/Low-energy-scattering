# -*- coding: utf-8 -*-
"""
Cria o esqueleto da entrada de HOJE no diário de desenvolvimento.

Uso (na raiz do projeto):
    python -m src.diario              # insere o template do dia no topo
    python -m src.diario --versao 1.6 # idem, já marcando a versão

Se já existir entrada com a data de hoje, não duplica (só avisa).
A versão atual é lida do pyproject.toml se não for passada.
"""
import argparse
import datetime
import os
import re

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIARIO = os.path.join(AQUI, "diario", "DIARIO.md")
PYPROJECT = os.path.join(AQUI, "pyproject.toml")

TEMPLATE = """## {data} — v{versao}

### Mudanças
- 

### Descobertas
- 

### Erros encontrados → lições
- 

### Estado no fim do dia
- 

### Próximos
- 

---

"""


def versao_atual():
    try:
        s = open(PYPROJECT, encoding="utf-8").read()
        m = re.search(r'version\s*=\s*"([^"]+)"', s)
        return m.group(1).rstrip(".0") if m else "?"
    except OSError:
        return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--versao", default=None, help="versão da entrada (ex.: 1.6)")
    args = ap.parse_args()

    hoje = datetime.date.today().isoformat()
    s = open(DIARIO, encoding="utf-8").read()
    if f"## {hoje}" in s:
        print(f"Já existe entrada de {hoje} — edite diario/DIARIO.md direto.")
        return

    v = args.versao or versao_atual()
    entrada = TEMPLATE.format(data=hoje, versao=v)
    # insere logo após o primeiro '---' (fim do cabeçalho)
    marca = "\n---\n"
    i = s.find(marca) + len(marca)
    s = s[:i] + "\n" + entrada + s[i:]
    open(DIARIO, "w", encoding="utf-8").write(s)
    print(f"Entrada de {hoje} (v{v}) criada em diario/DIARIO.md — preencha as seções.")


if __name__ == "__main__":
    main()

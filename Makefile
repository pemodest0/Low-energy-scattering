test:
	python -m pytest tests/ -q

figuras:
	python main.py --sem-ajustes

lit:
	python referencias/literatura.py

app:
	streamlit run app/Inicio.py

lock:
	pip freeze > requirements-lock.txt

.PHONY: test figuras lit app lock

sync:          ## sincroniza com a outra maquina (Mac <-> Windows)
	./infra/sincronizar.sh

mapa:          ## abre o mapa do projeto
	@cat MAPA.md

install:
	python -m pip install --upgrade pip &&\
		python -m pip install -r requirements.txt


setup: install


build-traits:
	python scripts/build_traits.py

parse-counts:
	python scripts/parse_counts.py

tune-assumptions:
	PYTHONPATH=src python scripts/tune_assumptions.py

evaluate-model:
	PYTHONPATH=src python scripts/evaluate_model.py


test:
	PYTHONPATH=src python -m pytest -vv --cov=growth_model tests/*


format:
	black src/*.py tests/*.py


lint:
	python -m pylint --disable=R,C src/growth_model.py


all: install test format lint
.PHONY: install setup build-traits parse-counts evaluate-model plot-backtest \
	holdout-ceiling results test format format-check lint all

install:
	python -m pip install --upgrade pip &&\
		python -m pip install -r requirements.txt


setup: install


build-traits:
	python scripts/build_traits.py

parse-counts:
	python scripts/parse_counts.py

evaluate-model:
	PYTHONPATH=src python scripts/evaluate_model.py

plot-backtest:
	python scripts/plot_backtest.py

holdout-ceiling:
	PYTHONPATH=src python scripts/holdout_ceiling.py

results: evaluate-model plot-backtest


test:
	PYTHONPATH=src python -m pytest -vv --cov=growth_model --cov=calibration tests/*


format:
	black src/*.py tests/*.py scripts/*.py

format-check:
	black --check src/*.py tests/*.py scripts/*.py


# Docstring checks (C0114-C0116) stay on; W0621 is pytest's fixture pattern.
lint:
	PYTHONPATH=src python -m pylint --disable=R,C --enable=C0114,C0115,C0116 src/*.py scripts/*.py
	PYTHONPATH=src python -m pylint --disable=R,C,W0621 --enable=C0114,C0115,C0116 tests/*.py


all: install format lint test

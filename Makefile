.PHONY: install setup build-traits fetch-rainfall parse-counts evaluate-model plot-backtest \
	holdout-ceiling results results-v2 results-v3 test format format-check lint all

install:
	python -m pip install --upgrade pip &&\
		python -m pip install -r requirements.txt


setup: install


build-traits:
	python scripts/build_traits.py

fetch-rainfall:
	python scripts/fetch_rainfall.py

parse-counts:
	python scripts/parse_counts.py

evaluate-model:
	PYTHONPATH=src python scripts/evaluate_model.py

plot-backtest:
	python scripts/plot_backtest.py

holdout-ceiling:
	PYTHONPATH=src python scripts/holdout_ceiling.py

results: evaluate-model plot-backtest

# Development species only; the confirmation run is a one-off (see README).
results-v2:
	PYTHONPATH=src python scripts/evaluate_state_space.py

results-v3:
	PYTHONPATH=src python scripts/survey_design.py
	python scripts/plot_survey_design.py
	PYTHONPATH=src python scripts/replay.py


test:
	PYTHONPATH=src:scripts python -m pytest -vv --cov=growth_model --cov=calibration --cov=state_space --cov=monitoring tests/*


format:
	black src/*.py tests/*.py scripts/*.py

format-check:
	black --check src/*.py tests/*.py scripts/*.py


# Docstring checks (C0114-C0116) stay on; W0621 is pytest's fixture pattern.
lint:
	PYTHONPATH=src python -m pylint --disable=R,C --enable=C0114,C0115,C0116 src/*.py scripts/*.py
	PYTHONPATH=src:scripts python -m pylint --disable=R,C,W0621 --enable=C0114,C0115,C0116 tests/*.py


all: install format lint test

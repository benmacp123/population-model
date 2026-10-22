.PHONY: install setup build-traits fetch-rainfall parse-counts evaluate-model plot-backtest \
	holdout-ceiling results results-v2 results-v3 lock requirements test format format-check \
	lint all


# uv manages the venv at venv/ (the portfolio's convention) rather than
# .venv; activate it before the other targets.
export UV_PROJECT_ENVIRONMENT = venv

# Sync venv/ to uv.lock exactly (all groups).
install:
	uv sync --frozen


setup: install


# Re-resolve uv.lock after editing pyproject.toml, then refresh the pip export.
lock:
	uv lock
	$(MAKE) requirements


# requirements.txt is generated for pip users; edit pyproject.toml instead.
requirements:
	uv export --frozen --no-emit-project --no-hashes --output-file requirements.txt > /dev/null


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

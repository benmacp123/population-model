install:
	python -m pip install --upgrade pip &&\
		python -m pip install -r requirements.txt


test:
	python -m pytest -vv --cov=growth_model tests/*


format:
	black *.py


lint:
	python -m pylint --disable=R,C growth_model.py


all: install test format lint
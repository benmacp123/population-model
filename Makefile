install:
	python -m pip install --upgrade pip &&\
		python -m pip install -r requirements.txt &&\
		python -m pip install -r requirements-dev.txt


test:
	python -m pytest -vv --cov=myproj test_hello.py


format:
	black *.py


lint:
	python -m pylint --disable=R,C growth_model.py


all: install test lint
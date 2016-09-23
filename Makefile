.PHONY: clean-pyc clean-build clean singularity

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "singularity - Creates singularity image"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -fr {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

singularity:
	@echo 'Creating Singularity image'
	mkdir -p dist
	sudo rm -f dist/*.img
	tar -c * > dist/txbr.tar
	vers=`cat VERSION | sed "s/\n//g"`; \
	echo 'version $vers'; \
	imgfile=`echo dist/txbr-$$vers.img` ; \
	whfile=`echo d3r-$$vers-py2.py3-none-any.whl` ; \
	echo 'image file $imgfile' ; \
	sudo singularity create -s 10000 $$imgfile ; \
	sudo singularity bootstrap $$imgfile singularity/txbrcentos.def $$vers; \
	rm -f dist/txbr.tar
	echo 'Singularity Image created $imgfile'
	ls -l dist

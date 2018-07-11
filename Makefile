.PHONY: clean-pyc clean-build clean singularity dist
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
        from urllib import pathname2url
except:
        from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
        match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
        if match:
                target, help = match.groups()
                print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -fr {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

singularity:  ## creates singularity image
	@echo 'Creating Singularity image'
	mkdir -p dist
	sudo rm -f dist/*.img
	tar -c * > dist/txbr.tar
	vers=`cat VERSION | sed "s/\n//g"`; \
	echo 'version $vers'; \
	imgfile=`echo dist/txbr.img` ; \
	whfile=`echo d3r-$$vers-py2.py3-none-any.whl` ; \
	echo 'image file $imgfile' ; \
	sudo singularity create -s 10000 $$imgfile ; \
	sudo singularity bootstrap $$imgfile singularity/txbrcentos.def $$vers; \
	rm -f dist/txbr.tar
	echo 'Singularity Image created $imgfile'
	ls -l dist

checkrepo: ## checks if remote repo is CRBS
	@therepo=`git remote get-url origin | sed "s/^.*://" | sed "s/\/.*//"` ;\
        if [ "$$therepo" != "nbcrrolls" ] ; then \
        echo "ERROR can only do a release from master nbcrrolls repo, not from $$therepo" ; \
        exit 1 ;\
        else \
        echo "Repo appears to be master nbcrrolls $$therepo" ; \
        fi

updateversion: ## Updates version by updating VERSION file
	@cv=`cat VERSION`; \
        read -p "Current ($$cv) enter new version: " vers; \
        echo "Updating VERSION with new version: $$vers"; \
        echo $$vers > VERSION

dist: clean ## creates distributable package
	@vers=`cat VERSION` ; \
        hvers=`cat VERSION | sed "s/\./-/g"` ;\
        txbrdirname=txbr-$$vers ;\
        distdir=dist/$$txbrdirname ;\
        /bin/mkdir -p $$distdir ;\
        cp -a * $$distdir/. ;\
        /bin/rm -rf $$distdir/dist ;\
	sed -i "s/txbr-stack-.*template/txbr-stack-$$hvers\&template/g" $$distdir/README.md ;\
        sed -i "s/dist\/txbr-v*\.img/dist\/txbr-$$vers\.img/g" $$distdir/README.md ;\
        sed -i "s/releases\/.*\/txbr.*\.json/releases\/$$vers\/txbr\_$$vers\_basic\_cloudformation.json/g" $$distdir/README.md ;\
	sed -i "s/txbr-source\/archive\/.*\.tar.gz/txbr-source\/archive\/v$$vers\.tar\.gz/"g $$distdir/README.md ;\
        sed -i "s/^tar -zxf txbr-source-.*tar.gz/tar -zxf txbr-source-$$vers.tar.gz/g" $$distdir/README.md ;\
        sed -i "s/^cd txbr-source-.*/cd txbr-source-$$vers/g" $$distdir/README.md ;\
        cat aws/basic_cloudformation.json | sed "s/@@VERSION@@/$${vers}/g" > dist/txbr_$${vers}_basic_cloudformation.json ;\
        tar -C dist/ -cz $$txbrdirname > $$distdir.tar.gz ;\
        ls -l dist

release: dist checkrepo ## package and upload a release to s3
	@echo "WARNING Creating new release in 15 seconds"
	@echo "If you dont want to do this hit Ctrl-c now!!!!"
	sleep 15
	@vers=`cat VERSION` ; \
        tarfile=txbr-$${vers}.tar.gz ;\
        cloudform=txbr_$${vers}_basic_cloudformation.json ;\
        aws s3 cp dist/$$cloudform s3://txbr-releases/$${vers}/$$cloudform --acl public-read ; \
        txbrdirname=txbr-$$vers ;\
        distdir=dist/$$txbrdirname ;\
        cp $$distdir/README.md . ;\
        branchy=`git branch --list | sed "s/^\* *//" | head -n 1` ;\
        git commit -m 'updated launch stack link' README.md ;\
        git push origin $$branchy ;\
        git tag -a v$${vers} -m 'new release' ; \
        git push origin v$${vers} ; \
        echo "Congratulations on the new release" ;\
        echo "A new tag v$${vers} has been pushed to github" ;\
        echo "Cloud formation file has been uploaded to s3://txbr-releases/$${vers}/$$cloudform"



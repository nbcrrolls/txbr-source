.. :changelog:

History
-------

3.1.2 (9-23-16)
---------------

* Switched code import of scipy.lib.blas to scipy.linalg.bias in
  txbr/util/mpfit/mpfit.py

* Fixed bugs in scripts/txbr_filter.py where wrong variable was used
  for path for .preali files and fixed typo in variable name used
  in line to set the path to the .remap file.

* Under examples/data fixed fhv6.tar.gz tarball by updating pixel sizing
  in fhv6b.preali cause it was set to 4 angstroms per pixel instead of
  1 (the value in fhv6a.preali)

* Updated setup.py to try to use setuptools first and if that fails to
  fallback on distutils.core

* Added Makefile and singularity target that creates a singularity txbr
  image.

* Added Vagrant configuration under vagrant/ directory

* Updated setup.cfg to be more generic.


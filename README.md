[imod]: http://bio3d.colorado.edu/imod/
[etomo]: http://bio3d.colorado.edu/imod/doc/etomoTutorial.html
[fhv6]: https://github.com/nbcrrolls/txbr-source/blob/master/examples/data/fhv6.tar.gz
[singularity]: http://singularity.lbl.gov/
[centos]: https://www.centos.org/
[python]: http://www.python.org/

# Transform based Tracking, Bundle adjustment and Reconstruction (TxBR)

Transform based Tracking, Bundle adjustment and Reconstruction (**TxBR**) is an electron tomography package developed on top of the [IMOD][imod] utilities. At this time, it does not offer any fiducial tracking options and it loosely follows the [ETomo][etomo] reconstruction scheme.

This source tree was derived from updated versions of **TxBR** code found here:

https://confluence.crbs.ucsd.edu/display/ncmir/TxBR

### Compatibility

* Works with [Python 2.6/2.7][python] on [Centos 6/7][centos]

### Dependencies

Simply put, a lot of software is needed to get this working. Below is a list of the software needed, but it is highly recommended that either the Virtual machine route or the [Singularity][singularity] route be used to use **TxBR** Both of these approaches are described in the [Installation](txbr-source#installation) section below.

* [IMOD][imod] should be installed in default location (/usr/local/IMOD)
* [opencv](http://opencv.org/) with development packages and python
* [Qt](https://www.qt.io) with development packages
* [Python 2.6/2.7][python] with development libraries
* [cln](http://www.ginac.de/CLN/) with development libraries
* [ginac](http://www.ginac.de) with utilities and development libraries
* [libtiff library](http://www.libtiff.org/)
* [numpy](http://www.numpy.org/)
* [scipy](https://www.scipy.org/)
* [python-matplotlib](http://matplotlib.org/)
* mpi4py-openmpi
* mpi4py-mpich
* [python-pillow](https://python-pillow.org/)
* canberra gtk3 library
* swig
* Cycler
* Pyrex 
* [swiginac](https://launchpad.net/ubuntu/+archive/primary/+files/swiginac_1.5.1.orig.tar.gz)
* gcc/g++

## Quickstart TxBR on the cloud

Click launch button to spin up the latest release of TxBR on the cloud (~20 minute spin up time):
**(Oregon region)**

[![Launch TxBR AWS CloudFormation link](https://s3.amazonaws.com/cloudformation-examples/cloudformation-launch-stack.png)](https://console.aws.amazon.com/cloudformation/home?region=us-west-2#/stacks/new?stackName=txbr-stack-3-1-2&templateURL=https://s3-us-west-2.amazonaws.com/txbr-releases/3.1.2/txbr_3.1.2_basic_cloudformation.json)

### Installation

[To simply try TxBR out follow these instructions](https://github.com/nbcrrolls/txbr-source/wiki/Using-TXBR-Vagrant-Virtual-machine)

**OR**
 
To build a [Singularity][singularity] image of TxBR do the following, assuming [Singularity][singularity] is installed:

```Bash
git clone https://github.com/nbcrrolls/txbr-source.git
cd txbr-source
make singularity
dist/txbr-v3.1.2-dev.img --help
```


**OR** 

Assuming all dependencies have been installed and **setup.cfg** paths are correct, then the following should work on Centos 7:

```Bash
git clone https://github.com/nbcrrolls/txbr-source.git
cd txbr-source
. /etc/profile.d/modules.sh
module load mpi/openmpi-x86_64
export LD_LIBRARY_PATH=/lib64:/usr/lib64/mpich/lib:/usr/local/IMOD/lib:/usr/local/IMOD/qtlib
python setup.py build
sudo python setup.py install
```

### Usage

As a starting point for reconstructing volumes, three files need to be provided for each series: preali files, rawtlt files and fid files (cf [IMOD][imod]). For instance, if we have two series called basenamea and basenameb, a total of six files should be accessible to the **TxBR** calculations: basenamea.preali, basenamea.fid, basenamea.rawtlt, basenamea.preali, basenamea.fid and basenamea.rawtlt.

It is very important to note that because **TxBR** jointly process multiple series, the fiducial files of all the disctinct series are related. The marker list should be the same in all the .fid files, even though one fiducial marker might be not visible in a certain series. 

To run the **TxBR** reconstruction, just type in from the command line: 

```Bash

runtxbr.py -b basename 

```

In the case of mutiple series of the same specimen, a comma separated list of the basenames has to be provided to the script runtxbr.py; for instance in case of a dual tilt series, one should type:

```Bash

runtxbr.py -b basenamea,basenameb 

```

Alignment of the micrographs, reorientation of the specimen reconstruction, filtering, backprojection and series combination will then follows in an automatic sequence; the final volume will then be displayed in an [IMOD][imod] windows. It is also possible to run manually each of the reconstruction steps. 

A series of options are available for the different **TxBR** scripts. Type:

```Bash

runtxbr.py --help

```

for a list of the options available to the runtxbr script. In particular, to start taking into account of the electron beam curvilinearity, you can set the n2 flag to an integer value higher than 1 (default value for a projective model); in practice, do not set a value larger than 6 for n2.

You can find a test dataset in the **examples/data** folder of the **TxBR** directory. Files needed to reconstruct a dual tilt series with **TxBR** are included in the tar file [fhv6.tar.gz][fhv6].

For more infomration visit the wiki:

https://github.com/nbcrrolls/txbr-source/wiki

# Bugs

Please report them [here](https://github.com/nbcrrolls/txbr-source/issues)


# License

See license in this file: [LICENSE.txt](LICENSE.txt)



# Credits & Acknowledgements

* Developers: Sebastien Phan, Alexander Ward Kulungowski, Raj Singh, Masako Terada, James Obayashi, Albert Lawrence

* This research benefitted from the use of credits from the National Institutes of Health (NIH) Cloud Credits Model Pilot, a component of the NIH Big Data to Knowledge (BD2K) program. 

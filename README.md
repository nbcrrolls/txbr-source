[imod]: http://bio3d.colorado.edu/imod/
[etomo]: http://bio3d.colorado.edu/imod/doc/etomoTutorial.html
[fhv6]: https://github.com/nbcrrolls/txbr-source/blob/master/examples/data/fhv6.tar.gz

# Transform based Tracking, Bundle adjustment and Reconstruction (TxBR)

Transform based Tracking, Bundle adjustment and Reconstruction (**TxBR**) is an electron tomography package developed on top of the [IMOD][imod] utilities. At this time, it does not offer any fiducial tracking options and it loosely follows the [ETomo][etomo] reconstruction scheme.

### Compatibility

* Works with Python 2.6, 2.7 on Linux environments

### Dependencies

Simply put, a lot of software is needed to get this working. 



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

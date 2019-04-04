<h1> Epylib  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;    <image src="https://raw.githubusercontent.com/mvilavidal/Epylib/master/images/signal.png" alt="signal" align="middle"  height="100" width="400"/> </h1>

[![DOI](https://zenodo.org/badge/179277414.svg)](https://zenodo.org/badge/latestdoi/179277414)

Epylib is an open-source Python 2 software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients.

## Description

The package is contained in the folder epylib. Basic documentation can be found in the wiki: https://github.com/mvilavidal/Epylib/wiki
.

Apart from the wiki, the script `tutorial.py` provides three comprehensive examples of the package functionality. For the first example, the script import simulated SEEG recordings from the folder `data/`. In the the second and third examples we analyze real power time courses during ictal epochs. The dataset and a detailed description can be found in OSF (https://osf.io/nmxfy/). Download the folder 'power_hdf5' and place it in the same folder that contains the file 'tutorial.py'

For a more detailed and techincal description of the algorithms implemented in this package, please refer to

Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. *Clinical Neurophysiology*, 128(6), 977-985. (https://doi.org/10.1016/j.clinph.2017.03.040)


## Dependencies

Epylib uses the HDF5 binary data format. Please install the h5py package (https://www.h5py.org/).

Additionally, epylib can automatically import recordings stored in EDF files. To enable this functionality the python library pyEDFlib needs to be installed (https://pyedflib.readthedocs.io).

## Installation

Download the folder epylib and add the containing directory in the search path (you can use the variable `sys.path`).


## Citing

If you use the source code, please make sure to reference both the package and the paper:

> Vila-Vidal, M. (2017). Epylib v1.0, https://github.com/mvilavidal/Epylib. Zenodo. (https://doi.org/10.5281/zenodo.2626639)

> Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. *Clinical Neurophysiology*, 128(6), 977-985. (https://doi.org/10.1016/j.clinph.2017.03.040)


## License

Epylib is a free open-source software released under the General Public License version 2.






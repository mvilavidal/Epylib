# Epylib <img src="https://raw.githubusercontent.com/mvilavidal/Epylib/master/images/signal.png" alt="signal" height="80" width="500" align="right" />

[![DOI](https://zenodo.org/badge/179277414.svg)](https://zenodo.org/badge/latestdoi/179277414) 


Epylib is an open-source Python software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients.

## Description

The package is contained in the folder epylib. Basic documentation can be found in the wiki: https://github.com/mvilavidal/Epylib/wiki. Besides, all functions and classes in the package contain Docstings information that can be accessed using the command `help()`.

The script `tutorial.py` provides three comprehensive examples of the package functionality. For the first example, the script import simulated SEEG recordings from the folder `data/`. In the the second and third examples we analyze real power time courses during ictal epochs. The dataset and a detailed description can be found in OSF (https://osf.io/nmxfy/). Download the folder 'power_hdf5' and place it in the same folder that contains the file 'tutorial.py'

For a more detailed and techincal description of the algorithms implemented in this package, please refer to:

> Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. *Clinical Neurophysiology*, 128(6), 977-985, https://doi.org/10.1016/j.clinph.2017.03.040.


## Dependencies

Epylib runs on Python 2. It uses the HDF5 binary data format to store array-like data with associated metadata. Please install the h5py package from https://www.h5py.org/.

Additionally, epylib can automatically import recordings stored in EDF files. To enable this functionality the python library pyEDFlib needs to be installed (https://pyedflib.readthedocs.io).

## Installation

Download the folder epylib and add the containing directory in the search path (you can use the variable `sys.path`).


## Citing

If you use the source code, please make sure to reference both the package and the paper:

> Vila-Vidal, M. (2019). Epylib v1.0, https://github.com/mvilavidal/Epylib. Zenodo, https://doi.org/10.5281/zenodo.2630604.

> Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. *Clinical Neurophysiology*, 128(6), 977-985, https://doi.org/10.1016/j.clinph.2017.03.040.


## License

Epylib is a free open-source software released under the General Public License version 2.






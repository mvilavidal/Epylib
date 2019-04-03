# Epylib

Epylib is an open-source Python 2 software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients.

## Dependencies

Epylib uses the HDF5 binary data format. Please install the h5py package (https://www.h5py.org/).

Additionally, epylib can automatically import recordings stored in EDF files. To enable this functionality the python library pyEDFlib needs to be installed (https://pyedflib.readthedocs.io).

## Description

The package is contained in the folder epylib. The script `tutorial.py` contains three examples showing the package functionality. For the first example, the script import simulated SEEG recordings from the folder `data/`. In the the second and third examples we analyze real power time courses during ictal epochs. The real data along with a description can be found in OSF (https://osf.io/nmxfy/). Download the folder 'power_hdf5' and place it in the same folder that contains the file 'tutorial.py'

Documentation is available online at x. For a more detailed description and applications, please refer to

Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. *Clinical Neurophysiology*, 128(6), 977-985. (https://doi.org/10.1016/j.clinph.2017.03.040)


## Installation

Download the folder epylib and add the containing directory in the search path (you can use the variable `sys.path`).

## Documentation

Documentation is available online at **.



## Citing

If you use the source code, please make sure to reference both the package and the paper:

x

Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. *Clinical Neurophysiology*, 128(6), 977-985. (https://doi.org/10.1016/j.clinph.2017.03.040)


## License

Epylib is a free open-source software released under the General Public License version 2.






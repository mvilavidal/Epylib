# -*- coding: utf-8 -*-
"""
Epylib is an open-source Python software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients.

Copyright (C) 2017, Manel Vila-Vidal
Contact details: m@vila-vidal.com / manel.vila-vidal@upf.edu

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License v2.0 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v2.0 for more details.

You should have received a copy of the GNU General Public License v2.0 along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""

#import numpy as np # 3-clause BSD license
#import sys
#import h5py # 3-clause BSD license
#import matplotlib.pyplot as plt # BSD-style license
#import scipy.io as sio # 3-clause BSD license
#from scipy import signal
#from scipy import fftpack
#try:
#    import pyedflib # 2-clause BSD license
#except ImportError:
#    print "Module pyEDFlib is not available. Loading data in EDF format is disabled."



from basics import *
from spectral import *
from methods import*
from stats import *
from plotting import *


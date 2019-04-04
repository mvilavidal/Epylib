# -*- coding: utf-8 -*-
"""
Epylib is an open-source Python software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients. Basic documentation can be found in the wiki: https://github.com/mvilavidal/Epylib/wiki.

If you use the source code, please make sure to reference both the package and the paper:
> Vila-Vidal, M. (2017). Epylib v1.0, https://github.com/mvilavidal/Epylib. Zenodo. (https://doi.org/10.5281/zenodo.2626639)

> Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. Clinical Neurophysiology, 128(6), 977-985. (https://doi.org/10.1016/j.clinph.2017.03.040)

-----------------------------------------------------------------------

Copyright (C) 2017, Manel Vila-Vidal
Contact details: m@vila-vidal.com / manel.vila-vidal@upf.edu

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License v2.0 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v2.0 for more details.

You should have received a copy of the GNU General Public License v2.0 along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""


import numpy as np


def med(a,axis=None):
    """
    Compute the median along the specified axis.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose median is desired.
    axis : int
        Axis or axes along which the means are computed. For details, see 'np.percentile'.
        If not specified the median of all values in 'a' is taken.
        
    Returns
    -------
    y : float or ndarray
        Returns the median or a new array containing the median values.
    
    """
    return np.percentile(a,50,axis=axis)


def nanmed(m,axis=None):
    """
    Compute the median along the specified axis, while ignoring nan values.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose median is desired.
    axis : int
        Axis or axes along which the means are computed. For details, see 'np.percentile'.
        If not specified the median of all values in 'a' is taken.
        
    Returns
    -------
    y : float or ndarray
        Returns the median or a new array containing the median values.
    
    """
    return np.nanpercentile(m,50,axis=axis)


def cohensd(x1,x2):
    """
    Compute the effect size of the mean difference between two datasets using the Cohen's d.
    
    Parameters
    ----------
    x1, x2 : ndarrays
        Arrays containing the rtwo datasets.
        
    Returns
    -------
    d : float
        Cohen's d.
    """
    n1=x1.size
    n2=x2.size
    s1=x1.std()
    s2=x2.std()
    s=np.sqrt(  (  (n1-1.)*s1**2  +  (n2-1.)*s2**2  )
                               / (n1+n2-2.)  )
    d=((x1.mean()-x2.mean())/s)
    return d


def norm(pw,m,s):
    """
    Z-score the values in pw with respect to the normal distribution N(m,s).
    
    Parameters
    ----------
    pw : array-like
        Array containing the values to normalise.
    m,s : float
        Values defining the mean and the standard deviation of the normal distribution N(m,s) with respect to which the values in pw are to be normalised.
        
    Returns
    -------
    pwN : arrray_like
        z-scored version of pw.
    """    
    pwN=np.zeros(pw.shape)
    pwN=(pw-m)/s
    return pwN



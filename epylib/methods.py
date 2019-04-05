# -*- coding: utf-8 -*-
"""
Epylib is an open-source Python software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients. Basic documentation can be found in the wiki: https://github.com/mvilavidal/Epylib/wiki.

If you use the source code, please make sure to reference both the package and the paper:
> Vila-Vidal, M. (2017). Epylib v1.0, https://github.com/mvilavidal/Epylib. Zenodo. (https://doi.org/10.5281/zenodo.2630604)

> Vila-Vidal, M., Principe, A., Ley, M., Deco, G., Campo, A. T., & Rocamora, R. (2017). Detection of recurrent activation patterns across focal seizures: Application to seizure onset zone identification. Clinical Neurophysiology, 128(6), 977-985. (https://doi.org/10.1016/j.clinph.2017.03.040)

-----------------------------------------------------------------------

Copyright (C) 2017, Manel Vila-Vidal
Contact details: m@vila-vidal.com / manel.vila-vidal@upf.edu

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License v2.0 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v2.0 for more details.

You should have received a copy of the GNU General Public License v2.0 along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""


import numpy as np # 3-clause BSD license
from basics import *
from stats import *
    

#### Artifact detector based on power and correlation between channels increase
    
def PxR_artifact_detector(pw,cuts=[60,-60],window=200,step=1):
    """
    Artifact detector using a sliding time window approach. The detector finds time windows where the product of the mean power and the correlation of all channels is two standard deviations above its median value in the selected period.

    Parameters
    ----------
    pw : MetaArray
        Array containing the instantaneous power of all channels. pw has shape ( number of channels , number of time points ).
    cuts : list, optional
        Time boundaries (in seconds) defining epochs where the detector is to be applied independently.
        Default is [60,-60].
    window : int, optional
        Window width (in time points) where the mean power and the correlation are comoputed.
        Default is 200.
    step : int, optional.
        Window step (in time points) of the sliding window approach.
        Default is 1.
        
    Returns
    -------
    timeskeep : ndarray
        1-D vector of 1 and 0. If timeskeep[i]==0, it means that an artifact affecting al channels occured at time stamp i and the signals are corrupted.

    """

    # Detect artifacts equally affecting all channels (detector function: correlation * mean power)
    R=[]
    P=[]
    j=0
    while(j+window<=pw.samples):
        aux=np.corrcoef(pw[:,j:j+window])
        indices=np.triu_indices_from(aux,k=1)
        aux=aux[indices]
        m=aux.mean()
        R.append(m)
        P.append(pw[:,j:j+window].mean())
        j+=step

    R=np.array(R)
    P=np.array(P)      
    Rc=R*P
    
    cuts=[cut*pw.sampling for cut in cuts]    
    for i in range(len(cuts)):
        if cuts[i]<0: cuts[i]=pw.samples+cuts[i]
    
    vlinfil=[False]*(window/2)
    von=window/2
    for bis in cuts:
        aux=Rc[von-window/2:bis-window/2]
        vlinfil=vlinfil+list(aux<np.percentile(aux,50)+2*aux.std())
        von=bis
    aux=Rc[von-window/2:]
    vlinfil=vlinfil+list(aux<np.percentile(aux,50)+2*aux.std())
    vlinfil+=[False]*(window/2-1)
    
    timeskeep=np.array(vlinfil,dtype='float')
    return timeskeep



def findNSE(nMA):
    """
    This function finds the normalized seizure ensemble among a set of seizures with nMA profiles specified in nMA.
    
    Parameters
    ----------
    nMA : ndarray or MetaArray
        Array of shape ( number of seizures, number of channels ), containing the normalized mean activation of each channel within each seizure.
    
    Returns
    -------
    outliersz : list
        The function returns 0 for a non-outlier seizure and 1 for outlier-seizures. Seizures for which most of the channels' nMA were NaN (probably due to artifacts) are assigned a -1 and can be considered as outliers too.      
    """
    
    m=nanmed(nMA,axis=0)
    s=np.nanstd(nMA,axis=0)
    # nans will be treated as outliers
    outliersz=[]
    for n in nMA:
        nanfrac=(n!=n).sum()/float(n.shape[0])
        if nanfrac>0.9:
            outliersz.append(-1)
        else:
            outlierchfrac=((n-m)>2*s).sum()/float(n.shape[0])
            if outlierchfrac>0.15:
                outliersz.append(1)
            else:
                outliersz.append(0)
    return outliersz



def findDS(nMA,usech=None):
    """
    This function finds the channel dense set (DS) from the nMA profiles that will be used to define a baseline distribution of inactivations.
    
    Parameters
    ----------
    nMA : ndarray or MetaArray
        Array of shape ( number of seizures, number of channels ), containing the normalized mean activation of each channel within each seizure.
        
    Returns
    -------
    DS : list
        List of channels contained in the dense set.
    """
    
    if usech is None: usech=range(nMAm.shape[0])
    u=usech
    Nch=len(u)
    
    grid=np.zeros((Nch,Nch-2))
    
    ffvar=np.inf
    for ch in u:
        fvar=np.inf
        O=u[:]
        O.remove(ch)
        I=[ch]
        k=0
        while(len(O)>0):
            v=np.inf
            for i in O:
                if(nMA[:,I+[i]].std()<v):
                    v=nMA[:,I+[i]].std()
                    tr=i
            O.remove(tr)
            I.append(tr)
            if len(I)<3: continue
            var=nMA[:,I].std()/len(I)
            grid[ch][k]=var
            k=k+1
            if var<fvar:
                fvar=var
                myI=I[:]
                fn=len(myI)
        if fvar<ffvar:
            ffvar=fvar
            ffn=fn
            fch=ch
            Isel=myI[:]    
    DS=Isel
    
    return DS





def getAO(pw,usech=None,uset=None, startt=60., endt=-60., activity_threshold=2, time_active=None, timerate_active=0.95):
    """
    This function computes the activation onset (AO) times for each channel from the (normalized)  power time courses in pw.
    
    Parameters
    ----------
    pw : MetaArray
        Array containing the instantaneous power of all channels. pw has shape ( number of channels , number of time points ).
    usech : list
        List of non-artifacted channels for which AO times are to be computed.
        By default all channels are considered.        
    uset : array-like of bools, optional
        If specified it must have length to the number of time points contained in pw. 'uset' is used as a mask, where a False identifies an artifacted time point. Artifacted time points will not be considered.
        By default all time points are considered as non-artifacted.
    startt : float, optional
        First time stamp (in seconds) to be considered for AO identification, typically the seizure onset.
        By default, 60.
    endt : float, optional
        Last time stamp (in seconds) to be considered for AO identification, typically the end of the seizure. Negative times indicate times counted from the end of the file into past.
        By deafult, -60.
    activity_threshold : float, optional
        Activation threshold above which a channel is considered to be active.
        Default is 2.
    time_active, timerate_active : floats, optional
        These parameters are jointly used to avoid the confounding effect of spurious short-lived activations and inactivations, respectively. A channel will be considered active when its activation crosses and remains above the specified threshold for at least a fraction 'timerate_active' (between 0 and 1) of the time points contained in a time window of length 'time_active' (in seconds). 
        Default values are 10 seconds and 0.95, respectively.
    
    Returns
    -------
    AO : ndarray
        Activation onset (AO) times for each channel.
    """
    
    N=pw.channels
    T=pw.samples
    fs=pw.sampling
    f=timerate_active
    
    B=pw>activity_threshold
    
    if usech is None: usech=range(N)
    if uset is None: uset=np.ones(T).astype(np.bool)
    n=int(10*fs) if time_active is None else int(time_active*fs)    

    t1=int(startt*fs)
    t2 = int(endt*fs)
    if t2<0: t2=T+t2
    
    AO = np.zeros(N)
    for ch in range(N):
        if ch in usech:
            for t in range(t1,t2):
                if uset[t] and ((B[ch,t:][uset[t:]][:n]).nonzero()[0].shape[0])>n*f:
                    AO[ch]=t
                    break
                AO[ch]=np.inf
        else:
            AO[ch]=np.inf
    AO=AO/float(fs)
    return AO

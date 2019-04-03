# -*- coding: utf-8 -*-
"""
Epylib is an open-source Python software for exploring, visualizing, and analyzing human intracranial EEG recordings from drug-resistant epilepsy patients.

Copyright (C) 2017, Manel Vila-Vidal
Contact details: m@vila-vidal.com / manel.vila-vidal@upf.edu

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License v2.0 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v2.0 for more details.

You should have received a copy of the GNU General Public License v2.0 along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""


import numpy as np # 3-clause BSD license
from scipy import signal, fftpack # 3-clause BSD license
from basics import *
    


def analytical(eeg):
    """
    Compute the analytical representation of a signal.

    Parameters
    ----------
    eeg : eegdata
        Eegdata object cobntaining the time series.     
        
    Returns
    -------
    xa : ndarray
        Analytical representation of eeg.data.
    """
    data=eeg.data
    n=eeg.samples   
    coeffs=fftpack.fft(data)
    coeffs[:,1:]*=2
    coeffs[:,int(np.ceil(n/2.)):]=0
    return fftpack.ifft(coeffs)



def analytical_bands(eeg,startf=None,endf=None,bandstep=.1,bandwidth=.1,adaptive=True):
    """
    Compute the analytical representations of each channel's eeg signal stored in 'eeg' in narrow-band frequency windows spanning a specified band-limited frequency range.

    Parameters
    ----------
    eeg : eegdata
        Eegdata object containing the time series.
    startf : float, optional
        Lower endpoint for the band-limited frequency range in which the analytical representations of the signal are to be computed.
        Default is 1.
    endf : float, optional
        Upper endpoint for the band-limited frequency range in which the analytical representations of the signal are to be computed.
        Default is fs/2.-fs/5., where fs stands for the sampling frequency.
    bandstep : float, optional
        Bandstep of the narrow-band frequency windows. The behavior of this parameter depends on the value of 'adaptive'.
        Default is 0.1.
    bandwidth : 
        Bandwidth of the narrow-band frequency windows. The behavior of this parameter depends on the value of 'adaptive'.
        Default is 0.1.
    adaptive : bool, optional
        This flag specifies whether the bandstep and bandwidth of the narrow bands should be kept constant or scaled proportional to the frequency of the band. If True, bandstep and bandwidth should be expressed as relative frequencies (relative to the lower bound of each narrow-band frequency range). If False, bandstep and bandwidth should be expressed as absolute frequencies.
        Default is True.
        
    Returns
    -------
    x_a_bands : ndarray
        Analytical representations of the signals stored in 'eeg' in the narrow bands defined by 'lbounds' and 'ubounds'. x_a_bands has shape (eeg.channels,Nf,eeg.samples), where Nf is the number of narrow-band frequency windows.
    lbounds : list
        Lower bounds of the narrow-band frequency windows spanning the specified frequency range.
    ubounds : list
        Upper bounds of the narrow-band frequency windows spanning the specified frequency range.
    """
    
    data=eeg.data
    n=eeg.samples
    sampling=eeg.sampling    

    if startf is None: startf=1.
    if endf is None: endf=sampling/2.-sampling/5.
  
    # if adaptive is True, step and bandwidth are relative to the frequency considered. For example:
    # if bandwidth=0.2, step=0.1, the partition will be 10-12,11-13.2,12.1-14.52,...,50.5-60.6,...
    lbounds=[]
    ubounds=[]
    l=startf
    if adaptive:
        u=(1+bandwidth)*l  
        while(u<endf):
            lbounds.append(l)
            ubounds.append(u)
            l=(1+bandstep)*l 
            u=(1+bandwidth)*l 
        if u<=sampling/2.:
            lbounds.append(l)
            ubounds.append(endf)
    else:
        u=l+bandwidth
        while(u<endf):
            lbounds.append(l)
            ubounds.append(u)
            l+=bandstep
            u=l+bandwidth
        if u<=sampling/2.:
            lbounds.append(l)
            ubounds.append(endf)

    bands=len(lbounds)
    coeffs=fftpack.fft(data)
    coeffs[:,1:]*=2  
    coeffs[:,int(np.ceil(n/2.)):]=0
    coeffsbands=np.tile(coeffs,(bands,1,1))
    coeffsbands=coeffsbands.swapaxes(0,1)
    freqs=fftpack.fftfreq(n,d=1./sampling)

    for i in range(bands):
        mask= (freqs>=lbounds[i])*(freqs<ubounds[i])
        coeffsbands[:,i,mask==False]=0
        
    x_a_bands=fftpack.ifft(coeffsbands)
   
    return x_a_bands, lbounds, ubounds



def spectrogram(eeg,ch=None,startf=None,endf=None,bandstep=.1,bandwidth=.1,adaptive=True):
    """
    Compute the spectrogram of the specified channel's signal using the Hilbert transform method. The signal is first filtered in narrow-band frequency windows spanning a specified band-limited frequency range. For each narrow band, the envelope is extracted using the Hilbert transform. Finally, the envelopes are squared to obtain instantaneous power time courses for each channel and narrow-band frequency window.

    Parameters
    ----------
    eeg : eegdata
        Eegdata object containing the time series.
    ch : int, optional
        Channel for which the spectrogram is to be computed. If no channel is specified, spectrograms for all channels' signals are computed at once.
    startf : float, optional
        Lower endpoint for the band-limited frequency range in which the analytical representations of the signal are to be computed.
        Default is 1.
    endf : float, optional
        Upper endpoint for the band-limited frequency range in which the analytical representations of the signal are to be computed.
        Default is fs/2.-fs/5., where fs stands for the sampling frequency.
    bandstep : float, optional
        Bandstep of the narrow-band frequency windows. The behavior of this parameter depends on the value of 'adaptive'.
        Default is 0.1.
    bandwidth : 
        Bandwidth of the narrow-band frequency windows. The behavior of this parameter depends on the value of 'adaptive'.
        Default is 0.1.
    adaptive : bool, optional
        This flag specifies whether the bandstep and bandwidth of the narrow bands should be kept constant or scaled proportional to the frequency of the band.
        If True, bandstep and bandwidth should be expressed as relative frequencies (relative to the lower bound of each narrow-band frequency range).
        If False, bandstep and bandwidth should be expressed as absolute frequencies.
        Default is True.
        
    Returns
    -------
    pw : MetaArray
        Spectrogram(s). If a channel is specified 'pw' has shape (Nf,eeg.samples), where Nf is the number of narrow-band frequency windows. If no channel is specified, 'pw' has shape (eeg.channels,Nf,eeg.samples). The metadata stored in eeg is passed to the MetaArray object 'pw', which has two additional attributes: pw.lbounds and pw.ubounds are lists containing the lower and upper bounds of the narrow-band frequency windows spanning the specified frequency range, respectively.
    """
    
    if ch is None:
        eeg_a,lbounds,ubounds=analytical_bands(eeg,startf=startf,endf=endf,bandstep=bandstep,bandwidth=bandwidth,adaptive=adaptive)
    else:
        eeg_a,lbounds,ubounds=analytical_bands(eeg.cut(channels=[ch]),startf=startf,endf=endf,bandstep=bandstep,bandwidth=bandwidth,adaptive=adaptive)
        eeg_a=eeg_a[0]
    pw=MetaArray(np.abs(eeg_a)**2)
    for key in eeg.__dict__:
        if key!='data':
            pw.__dict__[key]=eeg.__dict__[key]
    pw.lbounds=lbounds
    pw.ubounds=ubounds
    del eeg_a
    return pw
 

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
import matplotlib.pyplot as plt # BSD-style license
from basics import *


def tfplot(tf,startt=None,endt=None,vmin=None,vmax=None,perc=2,center=False, ticksevery=30, name=None,clabel='Power'):
    """
    Plot a spectrogram.
    
    Parameters
    ----------
    tf : MetaArray
        Spectrogram of shape ( number of frequency bands , time stamps ).
    startt : float, optional
        Initial time to show (in seconds).
        By default it is 0.
    endt : float, optional
        Last time to show (in seconds).
        By default it will be the last time point.
    vmin, vmax : float, optional
        These parameters are used to set the colorbar range. Values in 'tf' above and below 'vmin' and 'vmax' will be saturated to either value, respectively. If only 'vmin' is specified, 'vmax' is set to tf.max(). If only 'vmax' is specified 'vmin' is set to tf.min(). If neither parameter is specified (and center is False), they are automatically set to np.percentile(tf,perc) and np.percentile(tf,1-perc), respectively. However, if center is True, they are further adjusted so that the mean of 'tf' is centered.
    perc : bool, optional
        If either 'vmin' or 'vmax' is specified, this paramater is ignored. See the description of 'vmin' and 'vmax0 for details.
        Default is 2.
    center : bool, optional
        If either 'vmin' or 'vmax' is specified, this paramater is ignored. See the description of 'vmin' and 'vmax0 for details.
        Default is False.        
    ticksevery : float, optional
        Spacing between ticks in the time axis.
    name : string, optional
        If specified, the figure will be saved under the name 'name'. If not specified the figure is displayed but not saved.
    clabel : string, optional
        Set the label of the colorbar.
        Default is 'Power'.
    """
    
    if startt is None: startt=0
    if endt is None: endt=int(tf.span)
    if vmin is None and vmax is None:
        vmin=np.percentile(tf,perc)
        vmax=np.percentile(tf,100-perc)
        if center==True:
            aux=np.max((vmax-tf.mean(),tf.mean()-vmin))
            vmax=tf.mean()+aux
            vmin=tf.mean()-aux
    if type(tf[0,0])==np.bool_:
        plt.imshow(tf[:,int(startt*tf.sampling):int(endt*tf.sampling)],cmap='gray_r',interpolation='nearest',extent=(startt,endt,0,len(tf.lbounds)),origin='lower')
    else:
        plt.imshow(tf[:,int(startt*tf.sampling):int(endt*tf.sampling)],cmap='jet',vmin=vmin,vmax=vmax,interpolation=None,extent=(startt,endt,0,len(tf.lbounds)),origin='lower')
        
    if perc!=0: cbar=plt.colorbar(extend='both')
    else: cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=12)   
    cbar.set_label(clabel,size=14)
    
    plt.xticks(np.arange(startt,endt,ticksevery),fontsize=12)  
    plt.xlabel('time (s)',fontsize=14)
    plt.ylabel('frequency (Hz)',fontsize=14)
    ylabels=["{:3.3g}".format(tf.lbounds[i]) for i in range(len(tf.ubounds))]
    plt.yticks(np.arange(0,len(tf.lbounds),1)[::3],ylabels[::3],fontsize=12)
    
    plt.title('Spectrogram',fontsize=14)
    
    plt.axis('tight')
    if name is not None:
        plt.savefig(name,dpi=600)
        plt.close()
    else:
        plt.show() 
        

           
def tcplot(tc,usech=None,startt=None,endt=None,perc=2,name=None,clabel='Power',cmap='',ticksevery=60, usechlabels=False):
    """
    Colorplot showing the power carried by each channel as a function of time.
    
    Parameters
    ----------
    tc : MetaArray
        Array containing the power signal for each channel. 'tc' has shape ( number of channels , number of time points ).
    usech : list
        List of non-artifacted channels to use when defining the colorbar range.
        By default all channels are considered.
    startt : float, optional
        Initial time to show (in seconds).
        By default it is 0.
    endt : float, optional
        Last time to show (in seconds).
        By default it will be the last time point.
    perc : bool, optional
        The colorbar range is set from np.percentile(tc,perc) to np.percentile(tc,1-perc). Values in 'tc' out of this range will be saturated.
        Default is 2.
    name : string, optional
        If specified, the figure will be saved under the name 'name'. If not specified the figure is displayed but not saved.
    clabel : string, optional
        Set the label of the colorbar.
        Default is 'Power'.
    cmap : string, optional
        Colormap to use with plt.imshow.
        If not specified, the default colormap is used.       
    ticksevery : float, optional
        Spacing between ticks in the time axis.
    usechlabels : bool, optional
        If False, channels are labelled with an ordinal index. If True, labels stored in self.channellabels are used.
        Default is False.
    """
    
    fig,ax=plt.subplots(figsize=(10,6))
    if startt is None: startt=0
    if endt is None: endt=tc.span
    tc_w=tc[:,int(startt*tc.sampling):int(endt*tc.sampling)]      

    if usech is None: usech=range(tc.channels)
    vmin=np.percentile(tc_w[usech,:],perc)
    vmax=np.percentile(tc_w[usech,:],100-perc)
    if (type(tc[0,0])==np.bool_) or (cmap=='gray_r'):
        plt.imshow(tc_w, cmap='gray_r',interpolation='nearest',extent=(startt,endt,0,tc.channels),origin='lower')
    else:
        plt.imshow(tc_w,vmin=vmin,vmax=vmax,interpolation='nearest',extent=(startt,endt,0,tc.channels),origin='lower')

    if perc!=0: cbar=plt.colorbar(extend='both')
    else: cbar=plt.colorbar()
    cbar.set_label(clabel,size=14)
    cbar.ax.tick_params(labelsize=12)
    
    plt.xlabel('time (s)',fontsize=14)
    plt.xticks(np.arange(startt+ticksevery,endt,ticksevery),fontsize=12)
 
    plt.ylabel('channel',fontsize=14)
    if usechlabels:
        electrodes=[]
        chmark=[]
        for ch in range(tc_w.channels):
            electrode=''
            for c in tc_w.channellabels[ch]:
                if c.isalpha(): electrode+=c
            if electrode not in electrodes:
                electrodes.append(electrode)
                j=0
            if j in range(0,1000,5):
                chmark+=[ch]
            j+=1
        chmark=np.array(chmark)
        ylabels=[tc_w.channellabels[ch] for ch in chmark]
    else:
        ylabels=range(0,tc.channels,5)
        chmark=ylabels
    
    plt.yticks(chmark+0.5,ylabels,fontsize=12)

    plt.title('TCplot',fontsize=14)                

    plt.axis('tight')
    if name is not None:
        plt.savefig(name,dpi=600)
        plt.close()
    else:
        plt.show()


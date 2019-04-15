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


import epylib as epy
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt


# The script tutorial.py provides three comprehensive examples of the package functionality. For the first example, the script import simulated SEEG recordings from the folder 'data/'. In the the second and third examples we analyze real power time courses during ictal epochs. The dataset and a detailed description can be found in OSF (https://osf.io/nmxfy/). Download the folder 'power_hdf5' and place it in the same folder that contains the file 'tutorial.py'


#%% 1. IMPORT AND ANALYZE SIMULATED (S)EEG DATA FROM A FILE
### ---------------------------------------------

# This is simulated data


# Define patient and seizure ids
subject='P0'
seizure='sz0'
filename='P0sz0_example20signals.hdf5'


eeg=epy.eegdata()
eeg.loadHDF5('data/'+filename) # Import seeg recordings from a hdf5
print eeg.channels,eeg.samples
print eeg.sampling

eeg.notch_filter() # notch filter to remove the effect of AC

plt.figure(1)
eeg.plot([0,1,2,11,13],startt=10,endt=13, usechlabels=True) # plot some channels
plt.figure(2)
eeg.ps([1,2]) # power spectrum of channel 1


pw=epy.spectrogram(eeg,startf=1.,endf=150) # spectrogram between 1 and 150 Hz
epy.tfplot(10*np.log10(pw[5]),clabel='Power (dB)')

pwBroad=pw.sum(axis=1) # total power between 1 and 150 Hz
epy.tcplot(pwBroad,usechlabels=True) # tcplot

epy.saveHDF5(pwBroad,'data/'+'pw_broad_example20signals.hdf5')
pwBroad=epy.loadHDF5('data/'+'pw_broad_example20signals.hdf5')


#%% 2. ANALYZE ONE SEIZURE
### ----------------------




# The seizure epoch starts 60 s after the beginning of the recording and ends 60 s before the end of the recording.

## 2.1. EXAMPLE 1: BROADBAND
## -------------------------

pwBroad=epy.loadHDF5('power_hdf5/P1sz1_pwBroad.hdf5')
fs=pwBroad.sampling
N=pwBroad.channels
# Artifacted channels will not be considered. Here we define channels to be used henceforth:
u=range(N) # no channels are artifacted in P1 -> all channels will be used

# Detect artifacts equally affecting all channels (detector function: correlation * mean power)
timeskeep=epy.PxR_artifact_detector(pwBroad)
mask=timeskeep.astype(np.bool)

# Normalization of the pw_s with respect to a pre-ictal baseline obtained by pooling the power values of the first 40 seconds of all channels.
m=pwBroad[u,0:40*fs][:,mask[:40*fs]].mean()
s=pwBroad[u,0:40*fs][:,mask[:40*fs]].std()
pwBroadN=epy.norm(pwBroad,m,s)

# TCplot broadband
# ----------------
epy.tcplot(pwBroadN*mask,usech=u,clabel='z-score',usechlabels=True)

# MA, broadband, whole-seizure
# ----------------------------
plt.figure('MAbroadwhole')
MA=pwBroadN[u,60*fs:-60*fs][:,mask[60*fs:-60*fs]].mean(axis=1)
plt.plot(MA)
plt.title('Mean activation (MA)\n (broadband, whole seizure)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('MA')


# MA, broadband, 5 s
# ------------------
plt.figure('MAbroad5s')
MA=pwBroadN[u,60*fs:65*fs][:,mask[60*fs:65*fs]].mean(axis=1)
plt.plot(MA)
plt.title('Mean activation (MA)\n (broadband, 5s)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('MA')


## 2.2. EXAMPLE 2: HIGH-GAMMA
## --------------------------

pwHG=epy.loadHDF5('power_hdf5/P1sz1_pwHG.hdf5')
timeskeep=epy.PxR_artifact_detector(pwHG)
mask=timeskeep.astype(np.bool)
m=pwHG[u,0:40*fs][:,mask[:40*fs]].mean()
s=pwHG[u,0:40*fs][:,mask[:40*fs]].std()
pwHGN=epy.norm(pwHG,m,s)

# TCplot high-gamma
# -----------------        
epy.tcplot(pwHGN*mask,usech=u,startt=58,endt=62,clabel='z-score',usechlabels=True,ticksevery=1)

# MA, high-gamma, 100 ms
# ----------------------
plt.figure('MAhighgamma100ms')
MA=pwHGN[u,60*fs:int(60.1*fs)][:,mask[60*fs:int(60.1*fs)]].mean(axis=1)
plt.plot(MA)
plt.title('Mean activation (MA)\n (high-gamma, 0.1 ms)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('MA')




#%% 3. ACCUMULATE INFORMATION FROM DIFFERENT SEIZURES
### -------------------------------------------------


# Define clinically marked SOZ as the ground truth
soz=[20,21,22,32,33]

# Define the band to use: 'Broad', 'AB', 'DTh', 'G', 'HG'
bd='Broad'



## 3.1. AVERAGE NORMALIZED MEAN ACTIVATION (nMA) ACROSS SEIZURES
## -------------------------------------------------------------


# Find MA for all seizures:
MA=np.zeros((8,56))
Timeskeep=[] # will store the results of the artifact detector to use it later
for sz in range(8):
    pw=epy.loadHDF5('power_hdf5/P1sz'+str(sz+1)+'_pw'+bd+'.hdf5')
    fs=pw.sampling
    N=pw.channels
    u=range(N)
    timeskeep=epy.PxR_artifact_detector(pw)
    Timeskeep.append(timeskeep)
    mask=timeskeep.astype(np.bool)
    m=pw[u,0:40*fs][:,mask[:40*fs]].mean()
    s=pw[u,0:40*fs][:,mask[:40*fs]].std()
    pwN=epy.norm(pw,m,s)
    MA[sz,]=pwN[u,60*fs:-60*fs][:,mask[60*fs:-60*fs]].mean(axis=1)
    print 'seizure',sz,'done'
Timeskeep=np.array(Timeskeep)
MA=epy.MetaArray(MA,**pwN.__dict__)


# Plot all MA profiles together
plt.figure('MAsbroadwhole')
plt.plot(MA.transpose())
plt.title('Mean activation (MA)\n (broadband, whole seizure)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('MA')


# Normalize MA with respect to all channels' MA and find nMA:
m=MA[:,u].mean(axis=1)
s=MA[:,u].std(axis=1)
nMA=(MA-m[:,np.newaxis])/s[:,np.newaxis]
plt.figure('nMAsbroadwhole')
plt.plot(nMA.transpose())
plt.title('Normalized mean activation (nMA)\n (broadband, whole seizure)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('nMA')


# Find the NSE (normalized seizure ensemble)
outliersz=epy.findNSE(nMA[:,u])
NSE=[sz for sz in range(len(outliersz)) if outliersz[sz]==0]
print NSE


# Plot the average nMA over the NSE
nMAm=nMA[NSE,:][:,u].mean(axis=0)
nMAs=nMA[NSE,:][:,u].std(axis=0)
plt.figure('nMAavbroadwhole')
plt.errorbar(range(56),nMAm,yerr=nMAs,linewidth=2)
plt.title('Average nMA across the NSE\n (broadband, whole seizure)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('nMA')


# Compare average nMA inside and outside the SOZ. Use ranksums test and effect size (Cohen’s d)
print 'Clinically marked SOZ:', [str(MA.channellabels[ch]) for ch in soz]
nsoz=[ch for ch in u if ch not in soz]
pval=st.ranksums(nMAm[soz],nMAm[nsoz])[1]
effsize=epy.cohensd(nMAm[soz],nMAm[nsoz])
print 'Comparison of nMA between SOZ and non-SOZ sites yielded:'
print 'p-value:',pval
print "Effect size (Cohen's d):",effsize





## 3.2. AVERAGE ACTIVATION ONSET (AO) ACROSS SEIZURES
## --------------------------------------------------


# Find dense-set
DS=epy.findDS(nMA,usech=u)


# Plot the average nMA over the NSE and mark the DS in red
nMAm=nMA[NSE,:][:,u].mean(axis=0)
nMAs=nMA[NSE,:][:,u].std(axis=0)
plt.figure('DS')
plt.errorbar(range(56),nMAm,yerr=nMAs,linewidth=2,color='k')
plt.plot(DS,nMAm[DS],'ro')
plt.title('Average nMA across the NSE and dense set\n (broadband, whole seizure)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('nMA')


# Define power threshold in the NSE reference space
alpha=3 # z-score
m=nMA[:,DS].mean()
s=nMA[:,DS].std()
w_NSE=m+alpha*s


# Find AO for all seizures
AO=epy.MetaArray(np.zeros((8,56)),**MA.__dict__)
for sz in range(8):
    # Denormalize the threshold -> translate into each seizure's space
    M=MA[sz].mean()
    S=MA[sz].std()
    w=(w_NSE)*S+M # seizure-specific threshold
    
    pw=epy.loadHDF5('power_hdf5/P1sz'+str(sz+1)+'_pw'+bd+'.hdf5')
    fs=pw.sampling
    N=pw.channels
    u=range(N)
    timeskeep=Timeskeep[sz]
    mask=timeskeep.astype(np.bool)
    m=pw[u,0:40*fs][:,mask[:40*fs]].mean()
    s=pw[u,0:40*fs][:,mask[:40*fs]].std()
    pwN=epy.norm(pw,m,s)
    
    AO[sz]=epy.getAO(pwN,usech=u,uset=mask, startt=60., endt=-60., activity_threshold=w, time_active=2, timerate_active=0.95)
    print 'seizure',sz,'done'



# Plot all AO profiles together
AO[AO==np.inf]=np.nan
plt.figure('AOsbroad')
plt.plot(AO.transpose())
plt.title('Activation onset (AO)\n (broadband)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('AO')


# Plot the average AO over the NSE
AOm=np.nanmean(AO[NSE,:][:,u],axis=0)
AOs=np.nanstd(AO[NSE,:][:,u],axis=0)
plt.figure('AOavbroad')
plt.errorbar(range(56),AOm,yerr=AOs,linewidth=2)
plt.title('Average AO across the NSE\n (broadband)')
plt.xticks([0,10,20,32,44],[MA.channellabels[ch] for ch in [0,10,20,32,44]])
plt.xlabel('channel')
plt.ylabel('AO')


# Compare average AO inside and outside the SOZ. Use ranksums test and effect size (Cohen’s d)
print 'Clinically marked SOZ:', [str(MA.channellabels[ch]) for ch in soz]
nsoz=[ch for ch in u if ch not in soz]
pval=st.ranksums(AOm[soz],AOm[nsoz])[1]
effsize=epy.cohensd(AOm[soz],AOm[nsoz])
print 'Comparison of nMA between SOZ and non-SOZ sites yielded:'
print 'p-value:',pval
print "Effect size (Cohen's d):",effsize



## 3.3. PATIENT FEATURE MAP
## ------------------------


# Plot Patient map: (nMA,AO)
plt.figure('nMAvsAOwholebroad')
plt.scatter(nMAm,AOm)
plt.scatter(nMAm[soz],AOm[soz],marker='+',color='r')
for ch in range(N):
    plt.text(nMAm[ch],AOm[ch],AOm.channellabels[ch],fontsize=8)
plt.title('Average (nMA,AO) across the NSE\n (broadband)\nSOZ channels are marked with a red cross')
plt.xlabel('Average nMA')
plt.ylabel('Average AO')



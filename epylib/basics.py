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
import sys
import h5py # 3-clause BSD license
import matplotlib.pyplot as plt # BSD-style license
import scipy.io as sio # 3-clause BSD license
from scipy import fftpack, signal
try:
    import pyedflib # 2-clause BSD license
except ImportError:
    print "Module pyEDFlib is not available. Loading data in EDF format is disabled."
    


#### Class eegdata

class eegdata(object):
    """
    This is the basic class for (s)eeg data analysis.
    """   
    
    def  __init__(self,data=None,sampling=None,**kwargs):
        """
        Create an eegdata object.
        
        Parameters
        ----------
        data : array_like, optional
            A 2D array of shape (N,T) containing the timeseries of N recording sites along T time points.
        sampling : float, optional
            Sampling rate of the timeseries.
        **kwargs : optional
            *kwargs* are used to specify the initial set of parameters to be associated with the data. If dic is a dictionary containing the metadata it can be passed to the constructor using the double-star operator: eegdata(data, sampling, **dic).

            
        Returns
        -------
        eeg : eegdata
            The object containing the timeseries and the associated metadata.
        """
        
        if data is not None:
            assert (sampling is not None), "Specify sampling rate."
            self.data=np.array(data)
            self.channels=data.shape[0]
            self.samples=data.shape[1]
            self.sampling=sampling
            self.span=self.samples/float(self.sampling)
            self.__dict__.update(kwargs)
            
            
    def prescan_txtheader(self,source):
        """
        Read and print the header of a .txt file to check format.
        
        Parameters
        ----------
        source : string
            File name to be read
        """
        f=open(source, 'r')
        s=''
        for i in range(14):
            line=f.readline()
            if line[0:2]!='%%':
                s=s+line
        print s
    
    def __read_txtheader(self,f):
        line=f.readline()
        if line[0:2]!='%%': raise NameError('Check file format')
        line=f.readline()
        if line.split('\t',1)[0]!='% Archivo original:': raise NameError('Check file format')
        line=f.readline()
        if line.split('\t',1)[0]!='% Hora inicio/fin del archivo original:': raise NameError('Check file format')
        line=f.readline()
        if line.split('\t',1)[0]!='% Hora de inicio/fin del archivo exportado:': raise NameError('Check file format')        
        line=f.readline()
        if line.split('\t',1)[0]!='% Unidades:': raise NameError('Check file format')
        line=f.readline()
        if line[0:2]!='%%': raise NameError('Check file format')
        line=f.readline()
        if line.split('\t',1)[0]!='% Nombre del paciente:': raise NameError('Check file format')
        p=line.split()[4:]
        self.subject=p[1][0]+p[2][0]+p[0][0]
        line=f.readline()
        if line.split('\t',1)[0]!='% Tasa de muestreo': raise NameError('Check file format')
        self.sampling=int(float(line.split()[4]))
        line=f.readline()
        if line.split('\t',1)[0]!='% Canales': raise NameError('Check file format')
        self.totalchannels=int(line.split()[2])
        line=f.readline()
        if line.split('\t',1)[0]!='% ID del paciente ID del estudio': raise NameError('Check file format')
        line=f.readline()
        if line.split('\t',1)[0]!='% SN amplificador': raise NameError('Check file format')
        line=f.readline()
        if line[0]!='%': raise NameError('Check file format')
        line=f.readline()
        if line[0:9]!='% Format:': raise NameError('Check file format')
        line=f.readline()
        if line.split('\t',1)[0]!='% Fecha.Hora': raise NameError('Check file format')
        self.channellabels=np.array(line.split()[4:-1])
        self.channels=len(self.channellabels)
        line=f.readline()
        if line[0:2]!='%%': raise NameError('Check file format')
            
    def __read_txtdata(self,f,segmentlength_s=None):     
        data=[]
        eof=False
        i=0
        
        while(1):
            line=f.readline()
            i=i+1
            if line=='':
                eof=True
                break
            #if i==0:
            #    self.start=line.split()[0]+' '+line.split()[1]
            #data.append(line.split()[4:-1])
            data.append(line)
            #self.end=line.split()[0]+' '+line.split()[1]
            if segmentlength_s!=None:
                if i>=segmentlength_s*self.sampling:
                    break            
        if len(data)==0:
            return -1
        self.start=data[0].split()[0]+' '+data[0].split()[1]
        self.end=data[-1].split()[0]+' '+data[-1].split()[1]
        aux=np.genfromtxt(data,dtype=np.float32,invalid_raise=False)[:,4:-1]
        nanfrac=(np.count_nonzero(np.isnan(aux)))/float(aux.size)*100
        if nanfrac>=0.1:
            print 'WARNING:','{0:.1f}'.format(nanfrac),'% of numbers are nan and have been converted to 0.'
        self.data=np.nan_to_num(aux).transpose()
        self.samples=self.data.shape[1]
        self.span=self.samples/float(self.sampling)
        return eof
        
    def read_txt(self,source,start_s=0,read_s=None): 
        """
        Import data from .txt file. Designed for files with header in Spanish.

        Parameters
        ----------
        source : string
            File name to be read
        start_s : float (it must be a number such that start_s * samplig rate is an integer), optional
            Time span (in seconds) to be left out before starting reading.
            By deafault the method will start reading the first time sample contained in the file.
        read_s : float (it must be a number such that read_s * samplig rate is an integer), optional
            Time span (in seconds) to read.
            By default the method will read until the end of file is reached.
        """
        self.filename=source[:-4]     
        f=open(source, 'r')
        self.__read_txtheader(f)
        for i in range(start_s*self.sampling):
            f.readline()   
        eof=self.__read_txtdata(f,segmentlength_s=read_s)
        if eof==-1:
            print 'EOF was reached before reading could start.'
        elif eof==True:
            print 'EOF was reached. Time span read is: ', self.span, 's'
        else:
            print 'EOF not reached. Time span read is: ', self.span, 's'
        f.close()


    def txttomat(self,source,start_s=0,read_s=None,segmentlength_s=None):
        """
        Import data from .txt file and save in .mat. Designed for files with header in Spanish.
        This method allows to chunk the data and store in different .mat files.

        Parameters
        ----------
        source : string
            File name to be read
        start_s : float (it must be a number such that start_s * samplig rate is an integer), optional
            Time span (in seconds) to be left out before starting reading.
            By deafault the method will start reading the first time sample contained in the file.
        read_s : float (it must be a number such that read_s * samplig rate is an integer), optional
            Time span (in seconds) to read.
            By default the method will read until the end of file is reached.
        segmentlentgh_s: float (it must be a number such that start_s * samplig rate is an integer), optional
            Time span (in seconds) of each chunk.
            If not specified the whole recording will be stored in a single .mat file.
            
        """
        if segmentlength_s==None:
            self.read_txt(source,start_s=start_s,read_s=read_s)
            self.save_mat()
        else:
            self.filename=source[:-4]     
            f=open(source, 'r')
            self.__read_txtheader(f)
            for i in range(start_s*self.sampling):
                f.readline()
            
            filename=self.filename
            i=0
            totalread=0
            while(1):
                eof=self.__read_txtdata(f,segmentlength_s=segmentlength_s)
                if eof==-1:
                    break
                if self.span==0:
                    break
                self.filename=filename+'_'+str(i)
                self.save_mat()
                i=i+1
                totalread=totalread+self.span
                print 'file '+self.filename+'\t\t time span: '+str(self.span)
                if eof==True:
                    break
                if read_s!=None:
                    if totalread>=read_s:
                        break
                    elif read_s-totalread<segmentlength_s:
                        segmentlength_s=read_s-totalread
            if eof==False or eof==-1:
                print 'EOF not reached. Time span read: ', totalread, 's'
            else:
                print 'EOF was reached. Time span read: ', totalread, 's'


    def read_edf(self, filename, **kwargs):
        """
        Import data from .edf file.

        Parameters
        ----------
        filename : string
            File name to be read
        **kwargs : optional
            *kwargs* are used to specify the initial set of parameters to be associated with the data.
            If dic is a dictionary containing the metadata it can be passed to the constructor using the double-star operator: read_ef(filename,**dic).
        """
        if 'pyedflib' not in sys.modules: raise ImportError("Module pyEDFlib is not available. Loading data in EDF format is disabled.")
        f = pyedflib.EdfReader(filename)
        
        n=f.signals_in_file
        self.channels=n
        self.channellabels = f.getSignalLabels()
        N=f.getNSamples()[0]
        self.samples=N
        self.data = np.zeros((n, N))
        for i in np.arange(n):
            self.data[i, :] = f.readSignal(i)
        self.sampling=int(f.samplefrequency(0))
        self.span=N/float(self.sampling)
        self.filename=f.file_name
        self.start=f.getStartdatetime()
        
        self.__dict__.update(kwargs)
        
                
    def read_mat(self,filename):
        """
        Deprecated. Use loadMat instead.
        """
        self.loadMat(filename)

               
    def loadMat(self,filename):
        """
        Load dataset and associated metadata from .mat file.

        Parameters
        ----------
        filename : string
            File name to be read.
        """
        dic=sio.loadmat(filename,squeeze_me=True,mat_dtype=True)
        dic.pop('__globals__')
        dic.pop('__header__')
        dic.pop('__version__')
        for k in dic:
            setattr(self,k,dic[k])


    def save_mat(self,filename=None):
        """
        Deprecated. Use saveMat instead.
        """
        self.saveMat(filename)

    def saveMat(self,filename=None):
        """
        Save an eegdata object to file in .mat format.

        Parameters
        ----------
        filename : string, optional
            Output filename.
            If no filename is specified, the method will use self.filename.
        """
        if filename is None:
            filename=self.filename
        sio.savemat(filename,self.__dict__)
        
    
    def loadHDF5(self, filename=None):
        """
        Load dataset and associated metadata from hdf5 file.
        
        Parameters
        ----------
        filename : file-like object, string, optional
            Name of the file on disk, or file-like object.
            If no filename is specified, the method will use self.filename.     
        """
        
        if filename is None: filename=self.filename
        with h5py.File(filename, 'r') as hf:
            assert(len(hf.keys())==1), "The file contains more than one dataset."
            key=hf.keys()[0]
            dset=hf[key]
            self.data=dset[()]
            self.__dict__.update(dict(dset.attrs.items()))

        
    def saveHDF5(self,filename=None):
        """
        Save an eegdata object to file in .hdf5 format.
        
        Parameters
        ----------
        filename : file-like object, string, optional
            File or filename to which the data is saved.
            If no filename is specified, the method will use self.filename.
        """
        if filename is None: filename=self.filename
        if not filename.endswith('.hdf5'): filename+='.hdf5'
        with h5py.File(filename, 'w') as hf:
            dset=hf.create_dataset('eegdata',data=self.data)
            dic={}
            for key in self.__dict__:
                if key!='data': dic[key]=self.__dict__[key]
            dset.attrs.update(dic)
        
    

    def plot(self, channels=None, startt=None, endt=None, spacing = None, usechlabels=False, **kwargs):
        """
        Plot the (s)eeg signals.
        
        Parameters
        ----------
        channels : list, optional
            Selection of channels to display in the plot.
            By default all channels will be displayed. 
        startt : float, optional
            Initial time to show (in seconds).
            By default it is 0.
        endt : float, optional
            Last time to show (in seconds).
            By default it will be the last time of the signal.
        spacing : float, optional
            This parameter is used to control the spacing between signals.
        usechlabels : bool, optional
            If False, channels are labelled with an ordinal index. If True, labels stored in self.channellabels are used.
            Default is False.
        **kwargs : optional
            Optional parameters passed to plt.plot function
        """
        if channels is None: channels=range(self.channels)
        if startt is None: startt=0
        if endt is None: endt=self.span
        
        if 'linewidth' not in kwargs: kwargs['linewidth']=1.
        if 'color' not in kwargs: kwargs['color']='k'
        
        delta=0        
        for ch in channels:
            row=self.data[ch,int(startt*self.sampling):int(endt*self.sampling)]     
            delta+=-row.min() if spacing is None else 0
            plt.plot(np.arange(startt,endt,1./self.sampling), row + delta, **kwargs)
            string=str(self.channellabels[ch]) if usechlabels else str(ch)
            plt.text(startt+(endt-startt)*0.05,row[:int((endt-startt)*0.1*self.sampling)].max()+delta,string)            
            delta+=row.max()+1 if spacing is None else spacing            

        
    def cut(self, channels=None, startt=None, endt=None):
        """
        Extract a chunk of the signal for a subset of channels.
        
        Parameters
        ----------
        channels : list, optional
            Selection of channels.
            By default all channels will be taken. 
        startt : float, optional
            Initial time of the chunk (in seconds).
            By default it is 0.
        endt : float, optional
            Last time of the chunk (in secons).
            By default it will be the last time of the signal.
            
        Returns
        -------
        eeg : eegdata
            An eegdata object satisfying the specified requirements.
        
        """
        
        if channels is None:  channels=range(self.channels)             
        if startt is None:    startt=0
        if endt is None:  endt=self.span

        eeg=eegdata(self.data[channels,int(startt*self.sampling):int(endt*self.sampling)],
                                sampling=self.sampling,subject=self.subject)

        eeg.channellabels=self.channellabels[channels]
        #x.end='unk' if endt!=self.span else self.end
        #x.filename=self.filename+'*'
        #x.start='unk' if startt!=0 else self.start
        #x.totalchannels=self.totalchannels
        
        return eeg

        
    def notch_filter(self):
        """
        Apply a notch filter at 50 Hz and its multiples to remove the effect of AC artifacts.
        """
        self.data=notch(self.data,50.,self.sampling,rf=10.)
        for i in range(100,(self.sampling+2)/2,50):
            self.data=notch(self.data,float(i),self.sampling)
        if 'filter' not in self.__dict__: self.filter=''
        self.filter+='-notch'
        
    def bandpass(self,lowcut,highcut,order=2):
        """
        Apply a bandpass Butterworth filter forward and backwards to the signal.
    
        Parameters
        ----------
        lowcut : float
            Low cutoff frequency.
        highcut : float
            High cutoff frequency.            
        order : int, optional
            The order of the filter. As the filter is applied twice, the combined filter order will be twice that of the original.
            By deafult the order is 2.
        """
        self.data = butter_bandpass_filter(self.data, float(lowcut), float(highcut), float(self.sampling), order=order)
        if 'filter' not in self.__dict__: self.filter=''
        self.filter+='-bandpass('+str(lowcut)+','+str(highcut)+')'
    
        
    def ps(self,channels=None):
        """
        Plot the power spectrum of selected channels.
        
        Parameters
        ----------
        channels : int or list
            Selection of channels to display in the plot.
            By default all channels will be displayed.            
        """
        if type(channels)==int:
            channels=[channels]
        elif channels is None:
            channels=range(self.channels)
        for ch in channels:
            freqs,power=ps(self.data[ch],self.sampling)
            plt.plot(freqs,10*np.log10(power),label=str(ch))
            plt.xlabel('frequency (Hz)')
            plt.ylabel('power (dB)')
            plt.legend()
            plt.show()
            
    def downsample(self,newsampling):
        """
        Downsample the timeseries.
        
        Parameters
        ----------
        newsampling : int
            New sampling rate.    
        """
        if self.sampling/newsampling!=self.sampling/float(newsampling):
            print "New sampling rate must divide ",self.sampling
        else:
            rf=self.sampling/newsampling
            self.data=self.data[:,::rf]
            self.sampling=newsampling
            self.samples=self.data.shape[1]
            self.span=self.samples/float(self.sampling)


#### Class MetaArray (ndarray with associated dictionary) & basic I/O functions

class MetaArray(np.ndarray):
    """
    MetaArray is a subclass of the numpy ndarray class with extended metadata functionality.
    
    A MetaArray object has a set of associated attributes that store metadata. For a MetaArray marray, marray.__dict__ returns a dictionary with all attributes. The content of an attribute k can be accessed marray.k or marray.__dict__[k]
    
    Parameters
    ----------
    (for the __new__ method)
    
    input_array : array_like
        List or numpy array containing the data to be interpreted as a MetaArray type.
    **kwargs : optional
        *kwargs* are used to specify the initial set of parameters to be associated with the data. If dic is a dictionary containing the metadata it can be passed to the constructor using the double-star operator: MetaArray(input_array,**dic).
    """

    def __new__(cls, input_array, **kwargs):
        obj = np.asarray(input_array).view(cls)
        obj.__dict__=kwargs
        return obj
    

    def __array_finalize__(self, obj):
        # self is the new MetaArray
        # obj is the old MetaArray
        # For example, if we do b=a[:3], a is obj, while b is self
        if obj is None: return
        if getattr(obj, '__dict__', None) is None: return
        self.__dict__ = getattr(obj, '__dict__',None)

        
    def average(self,axis=None,weights=None,returned=False):
        """
        Compute the weighted average along the specified axis.
        
        Refer to `numpy.average` for full documentation.
        
        See Also
        --------
        numpy.average : equivalent function
        """
        
        out=myarray(np.average(self,axis=axis,weights=weights))
        out.__dict__=self.__dict__
        return out


# I/O: python dicionary <-> .mat
    
def loadMat(filename):
    """
    Load .mat file into a python dictionary.
    
    Parameters
    ----------
    filename : string
        Name of the file on disk.
        
    Returns
    -------
    dic: dictionary
        Dictionary containing the data imported from the .mat file.
    """ 
    
    dic=sio.loadmat(filename,squeeze_me=True,mat_dtype=True)
    dic.pop('__globals__')
    dic.pop('__header__')
    dic.pop('__version__')
    return dic

def saveMat(dic,filename):
    """
    Save a python dictionary to a file in .mat format.
    
    Parameters
    ----------
    dic: dictionary
        Python dictionary to be saved.
    filename : file-like object, string
        File or filename to which the data is saved.
    """
    sio.savemat(filename,dic)

# I/O: MetaArray <-> .npy file + .mat with metadata
    
def loadNpyMat(filename):
    """
    Load dataset from .npy file and associated metadata from .mat file.
    
    Parameters
    ----------
    filename : string
        Name of the file on disk without extension. The dataset and metadata will be read from the files filename+".npy" and filename+"_meta.mat", respectively.
        
    Returns
    -------
    marray: MetaArray
        MetaArray containing the imported array and the associated metadata.
    """ 
    marray=MetaArray(np.load(filename+'.npy'))
    dic=sio.loadmat(filename+'_meta.mat',squeeze_me=True,mat_dtype=True)
    dic.pop('__globals__')
    dic.pop('__header__')
    dic.pop('__version__')
    for k in dic:
        if k=='data':
            marray.__dict__['data_o']=dic[k]
        else:
            marray.__dict__[k]=dic[k]
    return marray

def saveNpyMat(marray,filename):
    """
    Save a MetaArray object to two separate files: a .npy file containing the dataset and a .mat file containing the metadata.
    
    Parameters
    ----------
    marray: MetaArray
        MetaArray to be saved.
    filename : string
        Name of the file on disk without extension. The dataset and metadata will be saved to filename+".npy" and filename+"_meta.mat", respectively.
    """

    dic=marray.__dict__
    np.save(filename+'.npy',np.array(marray))
    sio.savemat(filename+'_meta.mat',dic)

# I/O: MetaArray <-> .hdf5
       
def loadHDF5(filename):
    """
    Load dataset and associated metadata from hdf5 file.
    
    Parameters
    ----------
    filename : file-like object, string
        Name of the file on disk, or file-like object.
        
    Returns
    -------
    marray: MetaArray
        MetaArray containing the imported array and the associated metadata.
    """
    
    with h5py.File(filename, 'r') as hf:
        assert(len(hf.keys())==1), "The file contains more than one dataset."
        key=hf.keys()[0]
        dset=hf[key]
        marray=MetaArray(dset[()])
        marray.__dict__=dict(dset.attrs.items())
    return marray
        
def saveHDF5(marray,filename):
    """
    Save a MetaArray object to file in .hdf5 format.
    
    Parameters
    ----------
    marray: MetaArray
        MetaArray to be saved.
    filename : file-like object, string
        File or filename to which the data is saved.
    """
    
    if not filename.endswith('.hdf5'): filename+='.hdf5'
    with h5py.File(filename, 'w') as hf:
        dset=hf.create_dataset('metaarray',data=np.array(marray))
        dset.attrs.update(marray.__dict__)

# I/O: load metadata from a .hdf5 file to a python dictionary

def loadHDF5_metadata(filename):
    """
    Load metadata associated to a dataset stored in a hdf5 file.
    
    Parameters
    ----------
    filename : file-like object, string
        Name of the file on disk, or file-like object.
        
    Returns
    -------
    metadata: dict
        Dictionary containing the imported metadata.
    """

    with h5py.File(filename, 'r') as hf:
        assert(len(hf.keys())==1), "The file contains more than one dataset."
        key=hf.keys()[0]
        dset=hf[key]
        metadata=dict(dset.attrs.items())
    return metadata



#### Filters
    
def notch(data,interference_frequency,fs,rf=35.):
    """
    Apply a notch filter at a specified interference frequency.

    Parameters
    ----------
    data : array_like
        The array of data to be filtered.
    interference_frequency : float
        Interference frequency to be removed.     
    fs : float
        Sampling frequency of the data to which the designed filter is to be applied.      
    rf : float, optional
        Reducing factor of the notch filter.
        By deafult the reducing factor is 35.
        
    Returns
    -------
    y : ndarray
        Filtered data.
    """
    # w0 is the interference frequency expressed in cycles/half-cycle. Half-cycle corresponds to the Nyquist frequency
    # w0=1 for the Nyquist frequency (sampling/2)
    w0=interference_frequency/(fs/2.)
    b,a=signal.iirnotch(w0,rf)
    y=signal.lfilter(b,a,data)
    return y


def butter_bandpass(lowcut, highcut, fs, order=2):
    """
    Obtain the filter coefficients for a Butterworth bandpass digital filter.
    Refer to the function scipy.signal.butter for details.

    Parameters
    ----------
    lowcut : float
        Low cutoff frequency.
    highcut : float
        High cutoff frequency.        
    fs : float
        Sampling frequency of the data to which the designed filter is to be applied.      
    order : int, optional
        The order of the filter.
        By deafult the order is 2.
        
    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=2):
    """
    Apply a bandpass Butterworth filter forward and backwards to a signal.

    Parameters
    ----------
    data : array_like
        The array of data to be filtered.
    lowcut : float
        Low cutoff frequency.
    highcut : float
        High cutoff frequency.        
    fs : float
        Sampling frequency of data.      
    order : int, optional
        The order of the filter. As the filter is applied twice, the combined filter order will be twice that of the original.
        By deafult the order is 2.
        
    Returns
    -------
    y : ndarray
        The data filtered in the desired frequency band.
    """
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y
    

#### Other power-related functions
    
def ps(series, sampling):
    """
    Estimate power spectral density using the FFT.

    Parameters
    ----------
    series : array_like
        Time series.       
    sampling : float
        Sampling frequency of data.
        
    Returns
    -------
    f : ndarray
        Array of frequencies.
    power : ndarray
        Power spectrum of series. 
    """
    coeffs=fftpack.fft(series)#/series.shape[0]
    freqs=fftpack.fftfreq(series.shape[0],d=1./sampling)
    return abs(freqs[:freqs.shape[0]/2]), np.concatenate((abs(coeffs[0:1])**2,2*abs(coeffs[1:series.shape[0]/2-1])**2,abs(coeffs[series.shape[0]/2-1:series.shape[0]/2])**2))
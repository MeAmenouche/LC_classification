#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:14:34 2019

@author: melissa
"""
import numpy as np
import sncosmo
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

__all__ =['read_bandpass', 'plot_bands']

def read_bandpass(fname):
    """Read bandpass from two-column ASCII file containing wavelength and
    transmission in each line.
    """
    
    bands = sncosmo.get_bandpass(fname)
    d = {'Wavelength': bands.wave , 'Transmission': bands.trans} 
    
    data_filter = pd.DataFrame(data= d)
   
    return data_filter


def plot_bands(wave, trans, wave_R, trans_r, wave_g, trans_g, LabFilter):
    """Plot the bands each filter tabfilter with color and label(labFilter) and
        matplotlib object "ax" to be given, the figure will be created in the
        notebook
    """
   
    fig = plt.figure(figsize=[10,6])
    ax = fig.add_subplot(111)
    
    ax.plot(wave, trans,color = 'orange', label = LabFilter[0])
   
    
    ax.plot(wave_R, trans_r, color ='red', label = LabFilter[1])
    
    
    ax.plot(wave_g, trans_g, color = 'green', label = LabFilter[2])
    
    ax.set_xlabel('Wavelength')
    ax.set_xlabel('Filter transmission')
    ax.legend(loc = 'upper right')
    
    
    
    ax.set_xlim(3000, 10000)
    #ax.set_ylim(0, 0.001)
    #plt.yscale('log') 
    #line = ax.plot([0, 1100], [1, 1], color='black', linewidth=1, 
    #        linestyle="--")
    #ax.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.suptitle('Filter transmission bands', ha="center")
    plt.savefig('ZTf_filters')
    plt.show()
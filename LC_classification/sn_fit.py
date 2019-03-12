#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 17:10:19 2019
SNCOSMO FIT
@author: melissa
"""
import matplotlib.pyplot as plt
import numpy as np
#from ztfquery import marshal
#from ztfquery import io
import pandas as pd
import seaborn as nb
import sncosmo
from scipy import integrate



class Marshall_data:
    def __init__(self,ZTF_name):
        
        self.ZTF_name = ZTF_name
        
        
        

    def read(self):
        data = pd.read_csv("/Users/melissa/Data/ZTF/marshal/marshal.csv")
        
    def LC_data(self, ZTF_name):
        sourcefile = "/Users/melissa/Data/ZTF/marshal/lightcurves/"+source+"/marshal_lightcurve_"+source+".csv"
        lc = pd.read_csv(sourcefile)
        flag_mag = lc["magpsf"]<99 
        self.lc_mag = lc[flag_mag]
        flag_lim = lc["magpsf"]==99
        self.lc_lim = lc[flag_lim]

    def plot_lc(self, mag, lc_mag):
        flag_g = lc_mag["filter"]=='"g"'
        lc_g = self.lc_mag[flag_g]
        flag_r = self.lc_mag["filter"]=='"r"'
        lc_r = self.lc_mag[flag_r]
        flag_i = lc_mag["filter"]=='"i"'
        lc_i = lc_mag[flag_i]
        flag_glim = self.lc_lim["filter"] == '"g"'
        lc_glim = self.lc_lim[flag_glim]
        flag_rlim = self.lc_lim["filter"] == '"r"'
        lc_rlim = self.lc_lim[flag_rlim]
        flag_ilim = self.lc_lim["filter"] == '"i"'
        lc_ilim = self.lc_lim[flag_ilim]


        fig = plt.figure(figsize=[9,6])


        ax = fig.add_subplot(111)

        ax.errorbar(lc_g["jdobs"], lc_g["magpsf"], yerr=lc_g["sigmamagpsf"], 
                    marker='o', color="green", linestyle="None", ecolor="green"
                    , label='g band',)

        ax.errorbar(lc_r["jdobs"], lc_r["magpsf"], yerr=lc_r["sigmamagpsf"],
                    marker='o', color="red", linestyle="None", ecolor="red",
                    label='r band')

        ax.errorbar(lc_i["jdobs"], lc_i["magpsf"], yerr=lc_i["sigmamagpsf"], 
                    marker='o', color="orange", linestyle="None", 
                    ecolor="orange", label='i band')

        ax.legend(loc='best', shadow=True, fontsize='x-large')
        ax.scatter(lc_glim["jdobs"], lc_glim["limmag"], c="green", marker='+')
        ax.scatter(lc_rlim["jdobs"], lc_rlim["limmag"], c="red", marker='+')
        ax.scatter(lc_ilim["jdobs"], lc_ilim["limmag"], c="orange", marker='+')

        ax.set_xlabel("Time (day)")
        ax.set_ylabel("Magnitude") 
        plt.gca().invert_yaxis()

        
        
        
        
        
        
        
class SN_cosmo():
    
    #def __init__(self, )
    
    def register_ZTF_bands(ZTF_analysis_data = 'melissa/Data/ZTF/Filters/'):
        band_file_G = 'Filter_G.csv'
        band_file_I = 'Filter_I.csv'
        band_file_R = 'Filter_R.csv'
    
        filt2 = pd.read_csv(ZTF_analysis_data/band_file_G)
        wl_G = filt2['lambda']
        transmission_G = filt2['transmission']
        filt2 = pd.read_csv(ZTF_analysis_data/band_file_I)
        wl_I = filt2['lambda']
        transmission_I = filt2['transmission']
        filt2 = pd.read_csv(ZTF_analysis_data/band_file_R)
        wl_V = filt2['lambda']
        transmission_V = filt2['transmission']
    
        band_G = sncosmo.Bandpass(wl_G,transmission_G,wave_unit=u.AA,name='GZTF')    
        sncosmo.registry.register(band_G, force=True)
        band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='IZTF')    
        sncosmo.registry.register(band_I, force=True)
        band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RZTF')    
        sncosmo.registry.register(band_R, force=True)
        
        
    
    
def toto():
    print('Hello')

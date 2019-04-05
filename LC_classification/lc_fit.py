#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 16:14:31 2019

@author: melissa
"""
import numpy as np
import sncosmo
from .light_curve import MarshaLC
from .ZTF_bands import load_filters
from pathlib import Path
import pandas as pd

class Fit_lc(object):
    
    def __init__(self, path=Path('/Users/melissa/')):
        self.path = path
        load_filters()
        self.ztf_table = pd.read_csv(self.path/'Data/ZTF/marshal/marshal.csv')
        self.guy_sn = []
        self.rms_sn = []
        self.dic_dif = {} 
        self.t0_mod = []
    def read(self, sn_name):
        """
        TO DO
        """
        
        self.marshal = MarshaLC(sn_name, path=self.path)
        self.data = self.marshal.table_sncosmo()
        
        #self.ztf_table = pd.read_csv(self.path/'Data/ZTF/marshal/marshal.csv')
        flag = self.ztf_table['name'] == sn_name
        red = self.ztf_table[flag]
        s= red['redshift'][red.index]
        self.zhl  = float(s)
        r = red['ra'][red.index]
        self.ra = float(r)
        d = red['dec'][red.index]
        self.dec = float(d)
        self.mwebv = self.marshal.get_mwebv(self.ra, self.dec)
        #self.mwebv = self.marshal.mwebv
        
#    def crit_guy(self, sn_name, t_0):
#        """ To select SN with "good" light curves. We perform a first fit by setting
#        c to 0 to get the t_0. The fitted t_0 is an input to perform Guy selection
#        """
#        lc= self.marshal.curve
#        flag_mag = lc["magpsf"]<99 
#        lc_mag = lc[flag_mag] 
#        lc_mag['jdobs'] = lc_mag['jdobs'] - 2400000.5
#        
#        phase = (lc_mag['jdobs'] - t_0)/(1+self.zhl)
#        #print(phase)
#        cond1 = phase[(phase > -11) & (phase < 36)]  #Least 4 points
#        
#        cond2 = phase[(phase > -10) & (phase < 6) ]  #Least 1 point
#        
#        cond3 = phase[(phase > 6) & (phase < 21) ]   #Least 1 point
#        lc_r = lc_mag[(lc_mag['filter']== '"r"')]
#        lc_g = lc_mag[(lc_mag['filter'] == '"g"')]
#        phase_r = (lc_r['jdobs'] - t_0)/(1+self.zhl)
#        cond_r4 = phase_r[(phase_r > -8) & (phase_r < 10)]
#        phase_g = (lc_g['jdobs'] - t_0)/(1+self.zhl)
#        cond_g4 = phase_g[(phase_g > -8) & (phase_g < 10)]
#        
#        if (len(cond1) > 3 and len(cond2) > 0 and len(cond3) > 0 and 
#            len(cond_r4) >0 and len(cond_g4) > 0):
#                self.guy_sn.append(sn_name)
#                
        
            
#        return self.guy_sn
    def crit_guy(self, sn_name, t_0):
        """ To select SN with "good" light curves. We perform a first fit by setting
        c to 0 to get the t_0. The fitted t_0 is an input to perform Guy selection
        """
       
        self.read(sn_name)
        phase = (self.data['mjd'] - t_0)/(1+self.zhl)
        #print(phase)
        cond1 = phase[(phase > -11) & (phase < 36)]  #Least 4 points
        
        cond2 = phase[(phase > -10) & (phase < 6) ]  #Least 1 point
        
        cond3 = phase[(phase > 6) & (phase < 21) ]   #Least 1 point
        r = self.data['band']== "p48r"
        lc_r = self.data[r]
        g = self.data['band']== "p48g"
        lc_g = self.data[g]
        
        phase_r = (lc_r['mjd'] - t_0)/(1+self.zhl)
        
        cond_r4 = phase_r[(phase_r > -8) & (phase_r < 10)]
        phase_g = (lc_g['mjd'] - t_0)/(1+self.zhl)
        cond_g4 = phase_g[(phase_g > -8) & (phase_g < 10)]
        
        if (len(cond1) > 3 and len(cond2) > 0 and len(cond3) > 0 and 
            len(cond_r4) >0 and len(cond_g4) > 0):
                self.guy_sn.append(sn_name)
                
                
        return self.guy_sn
    
    def crit_rms(self, sn_name, t_0):
        """
        """
        mask = []
        for i in range(len(sn_name)):
            self.read(sn_name[i])
            print(i)
            r = self.data['band']== "p48r"
            lc_r = self.data[r]
            g = self.data['band']== "p48g"
            lc_g = self.data[g]
            phase_g = (lc_g['mjd'] - t_0[i])/(1+self.zhl)
            phase_r = (lc_r['mjd'] - t_0[i])/(1+self.zhl)
            cond_g = phase_g[(phase_g > -15) & (phase_g < 40)]
            cond_r4 = phase_r[(phase_r > -15) & (phase_r < 40)]
            max_r = lc_r['flux'] == np.max(lc_r['flux'])
            flu_r = lc_r[max_r]
            max_g = lc_g['flux'] == np.max(lc_g['flux'])
            flu_g = lc_g[max_g]
            
            if (len(flu_r) > 1 or len(flu_g) > 1):
                flu_r = flu_r[0]
                flu_g = flu_g[0]
            dif = abs(flu_r['mjd'] - flu_g['mjd'])
            
            if (len(cond_g)>3 and len(cond_r4)>3 and dif < 50):
            #if (len(cond_g)>3 and len(cond_r4)>3):
                self.dic_dif[sn_name[i]] = dif
                self.rms_sn.append(sn_name[i])
                mask.append(True)
            else:
                mask.append(False)
        self.mask = np.array(mask, dtype=bool)
        #np.ones°like(data,dtype=bool)
        return self.rms_sn, self.dic_dif, self.mask
        
        
        
    def rms_lc(self, sn_name):
        """ We calculate the rms of the flux of each SN to remove some SN not 
        removed by Guy selection such as the ones with plateau (no raise, no peak and no down) 
        """
        self.read(sn_name)
        r = self.data['band']== "p48r"
        lc_r = self.data[r]
        rms_r = np.sqrt(np.sum([x**2 for x in lc_r['flux']])/len(lc_r))
        
        g = self.data['band']== "p48g"
        lc_g = self.data[g]
        rms_g = np.sqrt(np.sum([x**2 for x in lc_g['flux']])/len(lc_g))
        
        return rms_r, rms_g
    
        
    def fitting(self, sn_name) :
        """Fit the light curve for a SN (data as Astropy table) using the model
        SALT2 z and mwebv are setted, it returns the values of the SALT2 parameters
        and their errors. We plot the fit on the notebook.
        """
        self.read(sn_name)
        self.source = sncosmo.get_source('salt2', version='2.4')
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, effects=[dust], 
                    effect_names=['mw'], effect_frames=['obs'])         
        self.model.set(mwebv=self.mwebv)
        self.model.set(z=self.zhl)
        #self.model.set(c=0)
        
        try:
            res, fitted_model = sncosmo.fit_lc(self.data, self.model, ['t0','x1',
                              'c', 'x0'], phase_range= (-15, 45))
            #res, fitted_model = sncosmo.fit_lc(self.data, self.model, ['t0','x1', 'x0'])
            #return res, fitted_model
        except:
        #except RuntimeError:
            print('Fail:',sn_name)
            res, fitted_model = 'fit fail', np.nan
            
        return res, fitted_model
    
    
    def fitting_z(self, sn_name):
        """SALT2 fit using SNCOSMO for the SN without => z free parameter of the
        fit
        """
        self.read(sn_name)
        self.source = sncosmo.get_source('salt2', version='2.4')
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, effects=[dust], 
                    effect_names=['mw'], effect_frames=['obs'])         
        self.model.set(mwebv=self.mwebv)
        #self.model.set(c=0)
        
        try:
            results, fit_model = sncosmo.fit_lc(self.data, self.model, ['t0','x1',
                              'c', 'x0'])
            #res, fitted_model = sncosmo.fit_lc(self.data, self.model, ['t0','x1', 'x0'])
            #return res, fitted_model
        except:
        #except RuntimeError:
            print('Fail:',sn_name)
            results, fit_model = 'fit fail', np.nan
            
        return results, fit_model
    
    
    
    
    def select_sn(self):
        """Extract only SN from the Marshall data
        """
        donnee = self.ztf_table
        
        flag_Ia = donnee['classification'] == 'SN Ia'
        flag_Ib = donnee['classification'] == 'SN Ib'
        flag_Ic = donnee['classification'] == 'SN Ic'
        flag_II = donnee['classification'] == 'SN II'
        flag_IIP = donnee['classification']== 'SN IIP'
        flag_91T = donnee['classification'] == 'SN Ia 91T-like'
        flag_91bg = donnee['classification'] == 'SN Ia-91bg'
        flag_IIn = donnee['classification'] == 'SN IIn'
        
        frame1 = donnee[flag_Ia]
        frame2 = donnee[flag_Ib]
        frame3 = donnee[flag_Ic]
        frame4 = donnee[flag_II]
        frame5 = donnee[flag_IIP]
        frame6 = donnee[flag_91T]
        frame7 = donnee[flag_91bg]
        frame8 = donnee[flag_IIn]
        frame = [frame1, frame2, frame3, frame4, frame5, frame6, frame7, frame8]
        self.sn_data = pd.concat(frame)
        self.sn_data = self.sn_data.drop(['Unnamed: 0', 'candid', 'iauname', 'release_auth', 'release_status'], axis=1)
        self.sn_data = self.sn_data.assign(t0=np.nan, t0_mod=np.nan, c=0, x0=np.nan, x1=np.nan,
                                           mwebv=np.nan, chi2=np.nan, ndof=np.nan, 
                                           chiperndof=np.nan, x1_err=np.nan, t0_err=np.nan,
                                           c_err=np.nan, x0_err=np.nan, x0_rel=np.nan,
                                           guy = np.nan)
        return self.sn_data
    
    
    
    def t_modif(self,sn_name, t0):
        """
        """
        for i in range(len(sn_name)):
            self.read(sn_name[i])
            r = self.data[self.data['flux'] != 0]
            #r_mag = self.data[r]
            p = t0[i] - r['mjd'][0]
            self.t0_mod.append(p)
        return self.t0_mod
        
   
    def write_fit(self, sn_name, to, x0, x1, c, mw_ebv, chi2, ndof,x1_err,to_err,
                  c_err, x0_err):
        """Write the parameters of the fits in the sn_data
        """
        
        flag_sn = self.sn_data['name'] == sn_name
        n = self.sn_data[flag_sn]
      
        self.sn_data.loc[n.index, 't0']  = to
        self.sn_data.loc[n.index, 'c'] = c
        self.sn_data.loc[n.index, 'x0'] = x0
        self.sn_data.loc[n.index, 'x1'] = x1
        self.sn_data.loc[n.index, 'mwebv'] = mw_ebv
        self.sn_data.loc[n.index, 'chi2'] = chi2
        self.sn_data.loc[n.index, 'x1_err'] = x1_err
        self.sn_data.loc[n.index, 't0_err'] = to_err
        self.sn_data.loc[n.index, 'ndof'] = ndof
        if ndof != 0:
            
            self.sn_data.loc[n.index, 'chiperndof'] = chi2/ndof
        
        self.sn_data.loc[n.index, 'c_err'] = c_err
        self.sn_data.loc[n.index, 'x0_err'] = x0_err
        self.sn_data.loc[n.index, 'x0_rel'] = x0_err/x0
        
        return self.sn_data
        
    
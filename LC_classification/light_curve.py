#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:02:34 2019

@author: melissa
"""

from astropy.table import Table
from astropy.io.ascii import InconsistentTableError
import numpy as np
import pandas as pd
from .ZTF_bands import _DEFAULT_FILTERS
from pathlib import Path
import sfdmap
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Galactic
import astropy.table.table
import warnings


class MarshaLC:
    """Class for the lightcurve of a single source in the Marshal
    Arguments:  
    name -- source name in the GROWTH marshal
    Options:
    ra          -- right ascension of source in deg
    dec         -- declination of source in deg
    sfd_dir     -- path to SFD dust maps if not set in $SFD_MAP
    user        -- Marshal username (overrides loading the name from file)
    passwd      -- Marshal password (overrides loading the name from file)
    filter_dict -- dictionary to assign the sncosmo bandpasses to combinations
                   of instrument and filter columns. This is only needed if there 
                   is non-P48 photometry. Keys are tuples of telescope+intrument 
                   and filter, values are the sncosmo bandpass names, 
                   see _DEFAULT_FILTERS for an example. 
    """
    def __init__(self, name, path=Path('/Users/melissa/') , ra=None, dec=None, redshift=None, classification=None,
                 mwebv=0., **kwargs):
       
        #kwargs = self._load_config_(**kwargs) #Vient de BaseTable à enlever
        self.path = path
        self.name = name
        self.redshift = redshift
        self.classification = classification
        self.filter_dict = kwargs.pop('filter_dict', _DEFAULT_FILTERS)
        
        
        if ra is not None and dec is not None:
            self.ra = ra
            self.dec = dec

        self.mwebv = mwebv
        
        # get the light curve into a table
        #path_way
        sourcefile = path/"Data/ZTF/marshal/lightcurves/"/name/("marshal_lightcurve_"+name+".csv")
        
        self.curve = pd.read_csv(sourcefile)
        self.clean_filter()
        dati = self.clean_lim()
        
        tab_pd = dati.as_matrix()
        
        
        #tab_pd = Table.from_pandas(curve)
        for i in range(len(tab_pd)):
            tab_pd[i,3] = tab_pd[i, 3].replace('"', '')
            tab_pd[i,1] = tab_pd[i, 1].replace('"', '')
            tab_pd[i,8] = tab_pd[i, 8].replace('"', '')
        
        #self.table_orig = Table.read(sourcefile, format='ascii.csv')
        self.table_orig = Table(rows = tab_pd, names =('Unnamed: 0', 'date',
                        'jdobs', 'filter', 'absmag', 'magpsf', 'sigmamagpsf',
                        'limmag', 'instrument'))
        
        
        
        
        self._remove_duplicates_()
        
        
    def clean_filter(self): #METTRE INSTRUMENT PAS FILTER
        self.light_curve = self.curve[(self.curve['instrument' ] == '"P48+ZTF"')]
        
                
        return self.light_curve
        
        
    def clean_lim(self):
        lc_r = self.light_curve.reset_index()
        lc_r = lc_r.drop(['index'],axis=1)
        r = lc_r['magpsf'] != 99
        r_mag = lc_r[r]
        indexes = r_mag.index.values.tolist()
        for i in range(indexes[0], indexes[len(r_mag)-1]):
            if lc_r['limmag'][i] != 99  and lc_r['magpsf'][i] == 99:
                lc_r = lc_r.drop(index = i)
                
                
            
            
        return lc_r
        
        
        
    
    
    def table_sncosmo(self):
        """Table of lightcurve data in the format sncosmo requires for fitting 
        """
        t = self.table
        
        
        zp = 25.0
        mag, magerr = t['magpsf'], t['sigmamagpsf']
        mjd   = t['jdobs'] - 2400000.5
        flux  = 10.**(-0.4*(mag-zp))
        eflux = flux * 0.4 * np.log(10.) * magerr
        
        #eflux = np.float64(eflux)
        zp = np.zeros(len(flux)) + zp
    
        mask = []
        zpsys = []
        band = []
        peakmag, peakmjd = 99,0.0
        for n,r in enumerate(t) :
            f = (r['instrument'], r['filter'])
            if r['magpsf'] > 90.: 
                flux[n] = 0.
                #eflux[n] = 10**(-0.4*(r['limmag']-zp[n]))/5.
                eflux[n] = np.float64(10**(-0.4*(r['limmag']-zp[n]))/5.)
                if eflux[n] == np.inf : 
                    eflux[n] = np.uint64(10**(-0.4*(r['limmag']-zp[n]))/5.)
                    #print(eflux[n])
                
                
                
            
            if f in self.filter_dict.keys():
                band.append(self.filter_dict[f])
                zpsys.append('ab')
                mask.append(True)
            else:
                    mask.append(False)
        
        mask = np.array(mask, dtype=bool)
        out = Table(data=[mjd[mask], band, flux[mask], eflux[mask], zp[mask], zpsys],
                names=['mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'])
        out.meta['z'] = self.redshift
        if self.mwebv is not None:
            out.meta['mwebv'] = self.mwebv
        
        return out    
    
    

    
        
        
    def _remove_duplicates_(self):
       """This function removes potential duplicates from the lightcurve,
       i.e. multiple entries for the same JD. If there are different values for the
       magnitude, the detections (mag < 99) are kept. Note that this may sometimes 
       leave two detections with different magnitudes at the same JD. 
       """
       t = self.table_orig
       mask = []
       t_obs = np.unique(t['jdobs'])
       for t_ in t_obs:
           if np.sum(t['jdobs'] == t_) == 1:
               mask.append(True)
           else:
               mags = t['magpsf'][t['jdobs'] == t_]
               if len(np.unique(mags)) == 1:
                    mask.append(True)
                    for k in range(len(mags) - 1):
                        mask.append(False)
               elif np.sum(np.unique(mags) < 90) == 1:
                   done = False
                   for m_ in mags:
                       if m_ < 90. and not done:
                           mask.append(True)
                           done = True
                       else:
                           mask.append(False)
               else:
                    mags_ = np.unique(mags)
                    mags_ = np.array(mags_[mags_ < 90])

                    done = [False for k in range(len(mags_))]
                    for m_ in mags:
                        if m_ < 90.:
                            k = np.where(mags_ == m_)[0][0]
                            if not done[k]:
                                mask.append(True)
                                done[k] = True
                        else:
                            mask.append(False)
       

       self.table = t[np.array(mask)]


    def get_mwebv(self, ra, dec):
        
        c = SkyCoord(ra, dec, unit= u.deg)
        m = sfdmap.SFDMap(self.path/'Data/sfddata-master', scaling=1.0)
        galcoords = c.transform_to(Galactic)
        return m.ebv(galcoords.l.deg, galcoords.b.deg)

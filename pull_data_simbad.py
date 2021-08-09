#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:48:25 2021

@author: ganesh
"""

import astropy
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
'''from astroquery.gaia import Gaia
width = u.Quantity(0.0001, u.deg)
height = u.Quantity(0.0001, u.deg)
coord = SkyCoord(ra=1.61057440603902, dec=64.1961580868088, unit=(u.degree, u.degree))

tb=Gaia.query_object(coord,width,height)
parallax=tb[0][9]
    spectral_type.append(result['SP_TYPE'][0])
    parallax.append(result['PLX_VALUE'][0])
    vsini.append(result['ROT_Vsini'][0])
    vmag.append(result['FLUX_V'][0])
'''
#for spectral type 
import pandas as pd
from astroquery.simbad import Simbad
#t=pd.read_table('/home/ganesh/Documents/project_ARIES/catlogue/IUCAA/total_list_of_downloaded_data.txt', header=None)

t=pd.read_table('/home/ganesh/Documents/project_ARIES/catlogue/IUCAA/CBe_names.csv', header=None)

spectral_type=[]
parallax=[]
parallax_error=[]
vmag=[]
bmag=[]
umag=[]
gmag=[]
jmag=[]
hmag=[]
kmag=[]
vmag_error=[]
bmag_error=[]
umag_error=[]
gmag_error=[]
jmag_error=[]
hmag_error=[]
kmag_error=[]
vsini=[]
vsini_error=[]
result=[]
name=[]
Teff=[]
logg=[]

bmagbib,umagbib,gmagbib,jmagbib,hmagbib, kmagbib=[],[],[],[],[],[]

spbib,plxbib,vmagbib,rotbib,teffbib=[],[],[],[],[]
c=np.array(t[0])

s = Simbad()
s.add_votable_fields('sptype')
s.add_votable_fields('rot')
s.add_votable_fields('fluxdata(V)')
s.add_votable_fields('parallax')
s.add_votable_fields('fe_h')
s.add_votable_fields('fluxdata(B)')
s.add_votable_fields('fluxdata(U)')
s.add_votable_fields('fluxdata(G)')
s.add_votable_fields('fluxdata(J)')
s.add_votable_fields('fluxdata(H)')
s.add_votable_fields('fluxdata(K)')

for i in range(len(c)):
    result.append(s.query_object(c[i]))
    name.append(c[i])
    spectral_type.append(result[i]['SP_TYPE'][0])
    spbib.append(result[i]['SP_BIBCODE'][0])
    Teff.append(result[i]['Fe_H_Teff'][0])
    logg.append(result[i]['Fe_H_log_g'][0])
    teffbib.append(result[i]['Fe_H_bibcode'][0])
    parallax.append(result[i]['PLX_VALUE'][0])
    parallax_error.append(result[i]['PLX_ERROR'][0])
    plxbib.append(result[i]['PLX_BIBCODE'][0])
    vsini.append(result[i]['ROT_Vsini'][0])
    vsini_error.append(result[i]['ROT_err'][0])
    rotbib.append(result[i]['ROT_bibcode'][0])
    vmag.append(result[i]['FLUX_V'][0])
    vmag_error.append(result[i]['FLUX_ERROR_V'][0])
    vmagbib.append(result[i]['FLUX_BIBCODE_V'][0])
    bmag.append(result[i]['FLUX_B'][0])
    bmag_error.append(result[i]['FLUX_ERROR_B'][0])
    bmagbib.append(result[i]['FLUX_BIBCODE_B'][0])
    
    umag.append(result[i]['FLUX_U'][0])
    umag_error.append(result[i]['FLUX_ERROR_U'][0])
    umagbib.append(result[i]['FLUX_BIBCODE_U'][0])
    gmag.append(result[i]['FLUX_G'][0])
    gmag_error.append(result[i]['FLUX_ERROR_G'][0])
    gmagbib.append(result[i]['FLUX_BIBCODE_G'][0])
    jmag.append(result[i]['FLUX_J'][0])
    jmag_error.append(result[i]['FLUX_ERROR_J'][0])
    jmagbib.append(result[i]['FLUX_BIBCODE_J'][0])
    hmag.append(result[i]['FLUX_H'][0])
    hmag_error.append(result[i]['FLUX_ERROR_H'][0])
    hmagbib.append(result[i]['FLUX_BIBCODE_H'][0])
    kmag.append(result[i]['FLUX_K'][0])
    kmag_error.append(result[i]['FLUX_ERROR_K'][0])
    kmagbib.append(result[i]['FLUX_BIBCODE_K'][0])
    
#    tb=np.vstack((name,spectral_type,parallax,vsini,vmag))
#    table=tb.transpose()
    table=np.column_stack((name,spectral_type,Teff,logg,parallax,parallax_error,
                           vsini,vsini_error,vmag,vmag_error,
                           bmag,bmag_error,umag,umag_error,gmag,gmag_error,
                           jmag,jmag_error,hmag,hmag_error,kmag,kmag_error,
                           spbib,plxbib,vmagbib,bmagbib,umagbib,gmagbib,
                           jmagbib,hmagbib,kmagbib,rotbib,teffbib))
    
    
    
np.savetxt('/home/ganesh/Documents/project_ARIES/catlogue/IUCAA/Bestar_Tess_spteffplxvsiniv.txt',table,
           newline= '\n',delimiter=',',
           header='Name,Spectral_type,Teff,logg,Parallax,Parallax_error,Vsini,vsini_error,Vmag,vmag_error,Bmag,bmag_error,Umag,umag_error,Gmag,gmag_error,Jmag,jmag_error,Hmag,hmag_error,Kmag,kmag_error,spbib,plxbib,vmagbib,bmagbib,umagbib,gmagbib,jmagbib,hmagbib,kmagbib,rotbib,teffbib' ,fmt='%s')

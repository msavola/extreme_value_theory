#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 08:35:08 2021

@author: mikkosavola
"""
#Script for plotting raw data from the Simms data and of the corresponding 
#dual distributions to see potential indication of Paretian distributions
#using the functions defined in file fat_tail_stats.py

import numpy as np
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import find_storms as fst
import time
import pickle
import fat_tail_stats as fts

#Set initial parameters
p_year=np.array([1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004],dtype=int) #,1996,1997,1998,1999,2000,2001
d_type='raw' #'raw' or 'dual'
big=9999
ms2_p=[1,2,3,4]
#Set maximum number of terms for summand in m2s plots
m2s_n_x1=1200
m2s_n_x2=1200
#Set k for Pickand's and Hill plots
k=2500
#Set cutoff point q for the smallest quantiles in the data to be discarded
q=0.0
#Set the relative random noise level
#(uniformly distributed) to be added to the data
noise=0


# def do_plots(x,a,ttl,zipf_xlbl,m2s_p):
#     ind=np.argwhere(a==p_year)
#     x_ind=x[ind]
#     print('\n')
#     print('Doing plots for'+ttl)
#     print('In year',str(p_year),'there are',str(len(x_ind)),'data points')
#     print('\n')
#     #Do me-plot
#     fts.me_plot(x,ttl)
#     #Do quantile-quantile plot
#     fts.qq_exp_plot(x,ttl)
#     #Do Zip plot
#     fts.zipf_plot(x,zipf_xlbl,ttl)
#     #Do maximum to sum plots
#     for ii in m2s_p:
#         fts.m2s_plot(x,ii,ttl)

#Read Sgeo and Sgr and eflux data
path="/Users/mikkosavola/Documents/Opiskelu_HY/phD/"\
                "data"
path.replace(" ", "")
path.replace("\n", "")
os.chdir(path)
file1=am.create_path(path,'ViRBO_ULF.txt')
file2=am.create_path(path,'e_flux.txt')
ulf,file1=am.reader_VIRBO(file1)
eflux,file2=am.reader_EFLUX(file2)

#Clean nans and -9999.9 values from the data
print('\nCleaning Sgeo')
sgeo,a_sgeo=fts.clean_data(ulf['Sgeo'],ulf['year'],big)
print('\nCleaning Sgr')
sgr,a_sgr=fts.clean_data(ulf['Sgr'],ulf['year'],big)
print('\nCleaning Fe1p2')
fe1p2,a_fe1p2=fts.clean_data(eflux['Fe1p2'],ulf['year'],big)
print('\nCleaning Fe130')
fe130,a_fe130=fts.clean_data(eflux['Fe130'],ulf['year'],big)
#l_year='2010'
#r_year=np.where(ulf['year']<l_year,1,0)

#Transform logged data into measurement values
sgeo=10**sgeo
sgr=10**sgr
fe1p2=10**fe1p2
fe130=10**fe130
#Find the q quantiles for the e fluxes
fe130_q=np.quantile(fe130,q)
fe1p2_q=np.quantile(fe1p2,q)   

#Print data of the run to the console
print('')
print('years'+str(p_year))
print('fe130 had originally '+str(len(fe130))+' data points')
print('fe1p2 had originally '+str(len(fe1p2))+' data points')
#Discard smallest qth tuantile of the e flux values and the corresponding
#year entries
a_fe130=a_fe130[fe130>fe130_q]
a_fe1p2=a_fe1p2[fe1p2>fe1p2_q]
fe130=fe130[fe130>fe130_q]
fe1p2=fe1p2[fe1p2>fe1p2_q]
print('fe130 has now '+str(len(fe130))+' data points')
print('fe1p2 has now '+str(len(fe1p2))+' data points')

#Add noise to avoid identical measurement values
if noise!=0:
    sgeo=fts.add_noise(sgeo,noise)
    sgr=fts.add_noise(sgr,noise)
    fe1p2=fts.add_noise(fe1p2,noise)
    fe130=fts.add_noise(fe130,noise)

#Remove duplicate values from Sgeo and Sgr
#dummy,sgeo_uniq_ind=np.unique(sgeo,return_index=True)
#dummy,sgr_uniq_ind=np.unique(sgr,return_index=True)
#sgeo=sgeo[sgeo_uniq_ind]
#sgr=sgr[sgr_uniq_ind]

#Do dual distributions
dual_sgeo=fts.dual_dist(sgeo)
dual_sgr=fts.dual_dist(sgr)
dual_fe1p2=fts.dual_dist(fe1p2)
dual_fe130=fts.dual_dist(fe130)

#Give plot names for figure files
if len(p_year)>1:
    year_0=p_year[0]
    year_1=p_year[-1]
    t_span=str(year_0)+'-'+str(year_1)
else:
    t_span=p_year
sgeo_ttl='sgeo_year_'+str(t_span)
sgr_ttl='sgr_year_'+str(t_span)
dual_sgeo_ttl='dual_sgeo_year_'+str(t_span)
dual_sgr_ttl='dual_sgr_year_'+str(t_span)
fe1p2_ttl='fe1p2_year_'+str(t_span)
dual_fe1p2_ttl='dual_fe1p2_year_'+str(t_span)
fe130_ttl='fe130_year_'+str(t_span)
dual_fe130_ttl='dual_fe130_year_'+str(t_span)

#Do plots for measurement data 
#The variable temp_label is used for setting x-label for Zipf plots
temp_lbl='Sgeo values (ULF spectral power)'
fts.do_plots(sgeo,k,sgeo_ttl,temp_lbl,ms2_p,a_sgeo,p_year,m2s_n_x1)
#temp_lbl='Dual of sgeo values (from log values)'
#fts.do_plots(dual_sgeo,dual_sgeo_ttl,temp_lbl,ms2_p,a_sgeo,p_year,m2s_n_x1)
temp_lbl='Sgr values (ULF spectral power)'
fts.do_plots(sgr,k,sgr_ttl,temp_lbl,ms2_p,a_sgr,p_year,m2s_n_x1)
#temp_lbl='Dual of sgr values (from log values)'
#fts.do_plots(dual_sgr,dual_sgr_ttl,temp_lbl,ms2_p,a_sgr,p_year,m2s_n_x1)
#Do plots for duals distributions
temp_lbl='Relativistic electron flux (keV−1cm−1s−1ster−1)'
fts.do_plots(fe1p2,k,fe1p2_ttl,temp_lbl,ms2_p,a_fe1p2,p_year,m2s_n_x2)
# temp_lbl='Dual of relativistic electron flux (from log values)'
# fts.do_plots(dual_fe1p2,dual_fe1p2_ttl,temp_lbl,ms2_p,a_fe1p2,p_year,m2s_n_x)
temp_lbl='Non-relativistic electron flux (keV−1cm−1s−1ster−1)'
fts.do_plots(fe130,k,fe130_ttl,temp_lbl,ms2_p,a_fe130,p_year,m2s_n_x1)
# temp_lbl='Dual of non-relativistic electron flux (from log values)'
# fts.do_plots(dual_fe130,dual_fe130_ttl,temp_lbl,ms2_p,a_fe130,p_year,m2s_n_x)





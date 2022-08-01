#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:28:00 2021

@author: mikkosavola
"""
import aux_mi as am
import main_mi as mm
import numpy as np
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import find_storms as fst
import time
import pickle
import fat_tail_stats as fts
import statsmodels.api as sm
import pandas as pd
#Write e flux data and the corresponding time stamps to file
#for a specified quantile

def flux_to_file(q,err):
    #Set initial parameters
    big=9999
   
    
    #Read Sgeo and Sgr and eflux data
    path="/Users/mikkosavola/Documents/Opiskelu_HY/PhD/EVT_flux_project/"\
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
    sgeo,d_sgeo=fts.clean_data(ulf['Sgeo'],ulf['day'],big)
    print('\nCleaning Sgr')
    sgr,a_sgr=fts.clean_data(ulf['Sgr'],ulf['year'],big)
    sgr.d_sgr=fts.clean_data(ulf['Sgr'],ulf['day'],big)
    print('\nCleaning Fe1p2')
    fe1p2,a_fe1p2=fts.clean_data(eflux['Fe1p2'],ulf['year'],big)
    fe1p2,d_fe1p2=fts.clean_data(eflux['Fe1p2'],ulf['day'],big)
    print('\nCleaning Fe130')
    fe130,a_fe130=fts.clean_data(eflux['Fe130'],ulf['year'],big)
    fe130,d_fe130=fts.clean_data(eflux['Fe130'],ulf['day'],big)
    #l_year='2010'
    #r_year=np.where(ulf['year']<l_year,1,0)
    
    #Transform logged data into measurement values
    sgeo=10**sgeo
    sgr=10**sgr
    fe1p2=10**fe1p2
    fe130=10**fe130
    
    #Set time stamps into a relative hourly timing, starting from 0
    for ii in range(len(a_fe130)):
        a_fe130[ii]=ii
        a_fe1p2[ii]=ii
        a_sgeo[ii]=ii
        a_sgr[ii]=ii

    #Find the q quantiles for the fluxes
    q=q
    fe130_qnt=np.quantile(fe130,q)
    fe1p2_qnt=np.quantile(fe1p2,q)  
    sgeo_qnt=np.quantile(sgeo,q)  
    sgr_qnt=np.quantile(sgr,q)  

    a_fe130=a_fe130[fe130>=fe130_qnt]
    a_fe1p2=a_fe1p2[fe1p2>=fe1p2_qnt]
    a_sgr=a_sgr[sgr>=sgr_qnt]
    a_sgeo=a_sgeo[sgeo>=sgeo_qnt]
    
    fe130_q=fe130[fe130>=fe130_qnt]
    fe1p2_q=fe1p2[fe1p2>=fe1p2_qnt]
    sgeo_q=sgeo[sgeo>=sgeo_qnt]
    sgr_q=sgr[sgr>=sgr_qnt]
    
    #Convert time stamps to seconds
    s_fe130=a_fe130#*60*60
    s_fe1p2=a_fe1p2#*60*60
    s_sgeo=a_sgeo#*60*60
    s_sgr=a_sgr#*60*60
    
    #Set errors, based on measurement precision, for Bayesian block analysis
    
    e_fe130=0.5*10**(err)*fe130_q #5*10^-5
    e_fe1p2=0.5*10**(err)*fe1p2_q #5*10^-5
    e_sgeo=0.5*10**(err)*np.ones(len(sgeo_q)) #5*10**(-3)
    e_sgr=0.5*10**(err)*np.ones(len(sgr_q)) #5*10**(-3)
    
    #Create the data frames
    f_fe130=pd.DataFrame({'TIME' : s_fe130,
                       'FLUX' : fe130_q,
                       'ERROR':e_fe130
                        })
    f_fe1p2=pd.DataFrame({'TIME' : s_fe1p2,
                       'FLUX' : fe1p2_q,
                       'ERROR':e_fe1p2                     
                        })
    f_sgeo=pd.DataFrame({'TIME' : s_sgeo,
                       'FLUX' : sgeo_q,
                       'ERROR':e_sgeo
                        })
    f_sgr=pd.DataFrame({'TIME' : s_sgr,
                       'FLUX' : sgr_q,
                       'ERROR':e_sgr
                        })
    
    
    #Write the data to files
    path=path+'/pkl'
    os.chdir(path)
    f_fe130.to_pickle(f"bb_input_log_fe130_q={q}_err={err}.pkl")
    f_fe1p2.to_pickle(f"bb_input_fe1p2_q={q}_err={err}.pkl")
    f_sgeo.to_pickle(f"bb_input_sgeo_q={q}_err={err}.pkl")
    f_sgr.to_pickle(f"bb_input_sgr_q={q}_err={err}.pkl")

err=-2.0
q=0.0
while q<0.1:
    while err<5:   
        flux_to_file(q,err)   
        err=err+0.5
    q=q+0.1
# q=0.95
# while q<0.999:
#     flux_to_file(q)
#     q=q+0.01
# q=0.995
# while q<0.9999:
#     flux_to_file(q)
#     q=q+0.001
# q=0.999
# while q<0.9999:
#     flux_to_file(q)
#     q=q+0.0001
    


    
    




    
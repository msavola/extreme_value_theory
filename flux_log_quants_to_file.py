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
import calc_errs as ce
#Write e flux data and the corresponding time stamps to file
#for a specified quantile


def flux_log_to_file(q,err,move_sd):
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
    
    #Set a running index for the flux values
    t_fe130=np.linspace(0,len(eflux['Fe130'])-1,len(eflux['Fe130']))
    t_fe1p2=np.linspace(0,len(eflux['Fe1p2'])-1,len(eflux['Fe1p2']))
    t_sgeo=np.linspace(0,len(ulf['Sgeo'])-1,len(ulf['Sgeo']))
    t_sgr=np.linspace(0,len(ulf['Sgr'])-1,len(ulf['Sgr']))
    
    
    #Clean nans and -9999.9 values from the data
    print('\nCleaning Sgeo')
    sgeo,a_sgeo=fts.clean_data(ulf['Sgeo'],ulf['year'],big)
    sgeo,d_sgeo=fts.clean_data(ulf['Sgeo'],eflux['doy'],big)
    sgeo,t_sgeo=fts.clean_data(ulf['Sgeo'],t_sgeo,1e10)
    print('\nCleaning Sgr')
    sgr,a_sgr=fts.clean_data(ulf['Sgr'],ulf['year'],big)
    sgr,d_sgr=fts.clean_data(ulf['Sgr'],eflux['doy'],big)
    sgr,t_sgr=fts.clean_data(ulf['Sgr'],t_sgr,1e10)
    print('\nCleaning Fe1p2')
    fe1p2,a_fe1p2=fts.clean_data(eflux['Fe1p2'],ulf['year'],big)
    fe1p2,d_fe1p2=fts.clean_data(eflux['Fe1p2'],eflux['doy'],big)
    fe1p2,t_fe1p2=fts.clean_data(eflux['Fe1p2'],t_fe1p2,1e10)
    print('\nCleaning Fe130')
    fe130,a_fe130=fts.clean_data(eflux['Fe130'],ulf['year'],big)
    fe130,d_fe130=fts.clean_data(eflux['Fe130'],eflux['doy'],big)
    fe130,t_fe130=fts.clean_data(eflux['Fe130'],t_fe130,1e10)
    #l_year='2010'
    #r_year=np.where(ulf['year']<l_year,1,0)
    
    # #Transform logged data into measurement values
    # sgeo=10**sgeo
    # sgr=10**sgr
    # fe1p2=10**fe1p2
    # fe130=10**fe130
    

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
    # s_fe130=a_fe130#*60*60
    # s_fe1p2=a_fe1p2#*60*60
    # s_sgeo=a_sgeo#*60*60
    # s_sgr=a_sgr#*60*60
    
    #Set errors, based on measurement s.d., for Bayesian block analysis   
    # e_fe130=0.5*10**(err)*fe130_q #5*10^-5
    # e_fe1p2=0.5*10**(err)*fe1p2_q #5*10^-5
    # e_sgeo=0.5*10**(err)*np.ones(len(sgeo_q)) #5*10**(-3)
    # e_sgr=0.5*10**(err)*np.ones(len(sgr_q)) #5*10**(-3)
    
    ns=0.5 #number of standard deviations to be used in calculating the error
    e_fe130=ce.calc_errs(fe130,err,ns,move_sd)
    e_fe1p2=ce.calc_errs(fe1p2,err,ns,move_sd)
    e_sgeo=ce.calc_errs(sgeo,err,ns,move_sd)
    e_sgr=ce.calc_errs(sgr,err,ns,move_sd)
    
    
    
    #Create the data frames
    f_fe130=pd.DataFrame({'TIME' : t_fe130,
                       'FLUX' : fe130_q,
                       'ERROR':e_fe130,
                       'YEAR':a_fe130 ,
                       'DOY':d_fe130
                        })
    f_fe1p2=pd.DataFrame({'TIME' : t_fe1p2,
                       'FLUX' : fe1p2_q,
                       'ERROR':e_fe1p2 ,
                       'YEAR':a_fe1p2 ,
                       'DOY':d_fe1p2
                        })
    f_sgeo=pd.DataFrame({'TIME' : t_sgeo,
                       'FLUX' : sgeo_q,
                       'ERROR':e_sgeo,
                       'YEAR':a_sgeo ,
                       'DOY':d_sgeo
                        })
    f_sgr=pd.DataFrame({'TIME' : t_sgr,
                       'FLUX' : sgr_q,
                       'ERROR':e_sgr,
                       'YEAR':a_sgr ,
                       'DOY':d_sgr
                        })
    
    if move_sd is not None:
        name_add="_move_sd"
    else:
        name_add=""
    
    #Write the data to files
    path=path+'/pkl'
    os.chdir(path)
    f_fe130.to_pickle(f"bb_input_log_fe130_err={err}_{ns}sigma_{name_add}.pkl")
    f_fe1p2.to_pickle(f"bb_input_log_fe1p2_err={err}_{ns}sigma_{name_add}.pkl")
    f_sgeo.to_pickle(f"bb_input_log_sgeo_err={err}_{ns}sigma_{name_add}.pkl")
    f_sgr.to_pickle(f"bb_input_log_sgr_err={err}_{ns}sigma_{name_add}.pkl")

err=20
move_sd=1
q=0.0
while q<0.1:
    while err<22:   
        flux_log_to_file(q,err,move_sd=1)   
        err=err+5
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
    


    
    




    
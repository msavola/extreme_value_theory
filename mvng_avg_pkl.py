#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 08:02:40 2022

@author: mikkosavola
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import pandas as pd
import os

#Moving averages of the data to pickle, data is x and averaging is over t
def mvng_avg_pkl(x,t,err,flux_name,move_sd=None):
    #Choose whether moving sd or non-moving sd data is used
    if move_sd is not None:
        name_add="_move_sd"
    else:
        name_add=""
    
    #Read the data
    f_name="bb_input_log_"+flux_name+f"err={err}{name_add}.pkl" 
    x=pd.read_pickle(f_name)
    
    #Calculate the moving averages
    mva=np.empty([0])
    for ii in range(t):
        DO JUST ONE AVERAGING WITH T BASED ON AUTOCORRLEATION
    
    
    #Save data to file
    mva.to_pickle(f"bb_input_log_{flux_name}_err={err}_mva_t={t}{name_add}.pkl")

path="/Users/mikkosavola/Documents/Opiskelu_HY/PhD/EVT_flux_project/"\
                    "data/pkl"
path.replace(" ", "")
path.replace("\n", "")
os.chdir(path)   
flux_name="fe130"   
start=1
stop=24*365 #one year in hours
lags=np.linspace(start,stop,num=start-stop+1)
errs=np.linspace(20,60,5,dtype=int)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:28:00 2021

@author: mikkosavola
"""
import aux_mi as am
import os
import pandas as pd
#Write the log of th e flux data and AL index to .pkl files


def flux_log_and_AL_to_file(): 
    
    #Read Sgeo and Sgr and eflux data
    path="/Users/mikkosavola/Documents/Opiskelu_HY/PhD/EVT_flux_project/"\
                    "data"
    path.replace(" ", "")
    path.replace("\n", "")
    os.chdir(path)
    file1=am.create_path(path,'OMNI_AL_1min_res_2013-2018.txt')
    file2=am.create_path(path,'e_flux.txt')
    AL,file1=am.reader_OMNI_AL(file1)
    eflux,file2=am.reader_EFLUX(file2)
    
    
    #Create the data frames
    f_fe130=pd.DataFrame({'HOUR' : eflux['hour'],
                       'FLUX' : eflux['Fe130'],
                       'YEAR':eflux['year'] ,
                       'DOY':eflux['doy']
                        })
    AL_ind=pd.DataFrame({'HOUR' :AL['hour'] ,
                       'AL' : AL['AL'],
                       'YEAR':AL['year'] ,
                       'DOY':AL['day'],
                       'MINUTE':AL['min']
                           })
    
    
    #Write the data to files
    path=path+'/pkl'
    os.chdir(path)
    f_fe130.to_pickle(f"log_fe130.pkl")
    AL_ind.to_pickle(f"AL_ind.pkl")



flux_log_and_AL_to_file()




    
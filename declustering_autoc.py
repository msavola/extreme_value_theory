#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 09:49:55 2022

@author: mikkosavola
"""

import os
import BayesBlocks as bb
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import gc
import calc_errs as ce

def declustering_autoc():
#Calculate the autocorrelation of the declustered data as the function of
#the declustering parameter cnt
    
    #Change file path
    path="/Users/mikkosavola/Documents/Opiskelu_HY/PhD/EVT_flux_project/"\
                    "data/pkl"
    path.replace(" ", "")
    path.replace("\n", "")
    os.chdir(path)

    #The input file is fixed because the error and sigma don't affect
    #the decluatering which uses only the flux values
    f_name="bb_input_log_fe130_err=20_4sigma__move_sd.pkl"
    df=pd.read_pickle(f_name)
    
    autos=np.empty([0])
    autos_ind=np.empty([0])
    h_max=24*100
    x_plot=np.linspace(1,h_max,h_max)
    
    u=3*10**5
    u_log10=np.log10(u)
    
    df2=df.loc[lambda df:df['FLUX']<7,:]
    flx2=np.array(df2.FLUX)
    flx=np.array(df.FLUX)
    
    
    autos2=np.empty([0])
    autos_ind2=np.empty([0])
    
    #Collect maxima and their indices
    maxis=np.empty([0])
    maxis2=np.empty([0])
    indis=np.empty([0])
    indis2=np.empty([0])
    
    
    
    #Cut flux and flux2 at the value >7
    flx_b=flx[:6615]
    flx2_b=flx2[:6615]
    
    flx_e=flx[6616:]
    flx2_e=flx2[6615:]
    
    #Pick the beginning or the end of flx and flx2
    #flx=flx_e
    #flx2=flx2_e
    
    
    
    for ii in range(1,h_max):
        #Find new maxima
        maxis,indis=ce.declustering(flx,u_log10,ii)
        maxis2,indis2=ce.declustering(flx2,u_log10,ii)
        #print(f"\nlen maxis={len(maxis)}")
        #print(f"len maxis2={len(maxis2)}\n")
        
        #if len(maxis2)!=len(maxis):
        #if ii==25:
        #      import pdb
        #      pdb.set_trace()
        
        #Plot the maxis and maxis2
        # plt.plot(indis,maxis,color='red',label='maxis')
        # plt.scatter(indis2,maxis2,color='blue',label='maxis2')
        # plt.legend()
        # plt.title('maxis and maxis 2')
        # plt.show()
        # maxis=np.concatenate(maxis,maxi*np.ones(1))
        # indis=np.concatenate(indis,indi*np.ones(1))
        # maxis2=np.concatenate(maxi2,maxi2*np.ones(1))
        # indis2=np.concatenate(indi2,indi2*np.ones(1))
        
        lm=len(maxis)
        #Break the loop if the length of maxis doesn't permit a further
        #iteration with a larger lag
        if lm<=ii:
            break
        #Convert log10 values to flux values
        maxis=10**maxis
        maxis2=10**maxis2
        if len(maxis)>0:
            #autoc=ce.my_autoc_lag(maxis,1)
            #autoc=autoc[0,1]
            #OR
            autoc=sm.tsa.acf(maxis,nlags=1)
            autoc=autoc[1]
            #OR autocovariance
            #autoc=sm.tsa.acovf(maxis,nlag=1)
            #autoc=autoc[1]
            
            autos_ind=np.concatenate([autos_ind,ii*np.ones([1])])
            autos=np.concatenate([autos,autoc*np.ones([1])])
            
            #autoc2=ce.my_autoc_lag(maxis2,1)
            #autoc2=autoc2[0,1]
            #OR
            autoc2=sm.tsa.acf(maxis2,nlags=1)
            autoc2=autoc2[1]
            #OR autocovariance
            #autoc2=sm.tsa.acovf(maxis2,nlag=1)
            #autoc2=autoc2[1]
            
            autos_ind2=np.concatenate([autos_ind2,ii*np.ones([1])])
            autos2=np.concatenate([autos2,autoc2*np.ones([1])])
            
            # if ii==40:
            #     import pdb
            #     pdb.set_trace()
        
    plt.plot(autos_ind,autos,label='autocorrelation',color='red')
    # g1=np.random.normal(loc=0,scale=1,size=120)
    # att=np.empty([0])
    # for ii in range(100):
    #     at=ce.my_autoc_lag(g1,ii)
    #     at=at[0,1]
    #     att=np.concatenate((att,np.ones(1)*at))
    # plt.plot(np.linspace(2,100,99),att[1:])
    # plt.show()
        
    plt.plot(autos_ind2,autos2,label='autocorrelation (log flx<7)',color='blue')
    const=0.1
    const_vec=const**np.ones(len(x_plot))
    plt.plot(x_plot,const_vec,label=f"autocorrelation={const}",color='black')
    plt.title("Autocorrelation with lag 1")
    plt.ylabel("autocorrelation")
    plt.xlabel("declustering parameter (h)")
    plt.xlim([0,lm])
    plt.legend()
    plt.savefig(f"my_autoc_fe130_u={u}_declustering.pdf")
    
declustering_autoc()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:28:00 2021

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
import scipy.stats.distributions as ssd
import fat_tail_stats as fts

#Fit a generalised Pareto distribution to the given data
#and get the distribution's parameters as a result
#Declustering is done using Bayesian blocks (J. Scargle)



            
    

#Calculate a moving average of vector x over time span d
def move_avg(flux,d):
    #Convert d from hours to days
    d=int(d*24)
    
    avgs=np.empty([0])
    lf=len(flux)
    
    #Check that d<lf
    if d>=lf:
        print("d must be less than the length of the flux")
        return
    
    #Calculate the first d/2 data points    
    d_flr=int(np.floor(d/2))
    subs=flux[:d_flr]
    avg=np.mean(subs)
    avgs=np.concatenate((avgs,np.ones([len(subs)])*avg)) 
    
    
    #Do averages for the data points from [:d/2] onwards, taking
    #the nearest d values around the data point
    for ii in range(lf-d):
        subs=flux[ii:ii+d]
        avg=np.mean(subs)
        avgs=np.concatenate((avgs,np.ones([1])*avg))
    #Do the averages for the final d/2 data points
    subs=flux[lf-d_flr:]
    avg=np.mean(subs)
    avgs=np.concatenate((avgs,np.ones([len(subs)])*avg)) 
    
    return avgs
    

#def plot_GPD_bb(q,errs,flux_name,plt_GPD=None,slice_df=None,t_add=None,move_sd=None,cut_first=None):
def plot_GPD_bb(errs,flux_name,plt_GPD=None,slice_df=None,t_add=None,move_sd=None,cut_first=None,mv_avg=None,fp_rate=0.003,thr=None):
    #Loop indices
    jj=0
    kk=0
    #Number of standard deviations to be usefld in calculating the error
    ns=4
    #Quantile for error calculation
    iq=0.25
    #Threshold value for GPD fitting based on me plots
    u=3*10**5
    log_u=np.log10(u) #set according to the me plot (use log because flux is log)
    
    #Use cut_first to decide whether the thresholding is done before the BB analysis
    if cut_first is not None:
        cut_txt="_cut_first"
    else:
        cut_txt=""
    #Choose whether a moving average is used or not
    if mv_avg is not None:
        avg_txt="_move_avg="
        #Error 42 is the about the smallest"zero autocorrelation time"
        avg_d=42
    else:
        avg_txt=""
        avg_d=""


    
    #Figure for plotting the BB data
    #######
    fig12=plt.figure(12)
    #Create an nxn figure with subplots based on the length of errs
    n_subpl=np.ceil(np.sqrt(len(errs)))
    n_subpl=int(n_subpl)
    fig12,ax12=plt.subplots(n_subpl,n_subpl,sharey='all',sharex='all',squeeze=False)
    #######
    
    if move_sd is not None:
        name_add="_move_sd"
    else:
        name_add=""
    
    for err in errs:
        
        #######
        f_name="bb_input_log_"+flux_name+f"_err={err}_{ns}sigma_{name_add}.pkl"
        df=pd.read_pickle(f_name)
        
        if slice_df is not None:
        #Slice the datafrane
            df=df.loc[lambda df:df['YEAR']==2003,:]
            df=df.loc[lambda df:df['DOY']>=244,:] #300, 1, 324, Oct 1st=274
            df=df.loc[lambda df:df['DOY']<=400,:] #306, 59, 324
            
        if t_add is None:
            #Title add-on, a date, for specific storms
            t_add=""
        
        #Do thresholding before the BB analysis, and set the values below the 
        #threshold to the threshold
        #u=10**np.quantile(df.FLUX,q)
        if cut_first is not None:
            df.FLUX[df.FLUX<log_u]=log_u
            df.ERROR=ce.calc_errs(df.FLUX,err,ns,move_sd)
        
        #Use moving averages, if desired
        if mv_avg is not None:
            df.FLUX=move_avg(df.FLUX,int(avg_d))
            #Also recalculate the errors
            df.ERROR=ce.calc_errs(df.FLUX,avg_d,ns,move_sd)
        
        #Remove the suspicious flux value of the order 10**7 in the fe130 data
        if flux_name=='fe130':
            df=df.loc[lambda df:df['FLUX']<7,:]
            t_add=t_add+"flux<7"
        
        #Do declustering of the data
        thr=thr
        log_maxis,indis=ce.declustering(np.array(df.FLUX[:]),log_u,thr)
        indis=indis.astype(int)
        #Convert maxis from log10 of flux to flux
        maxis=10**log_maxis
        
        #####UNCOMMENT STARTING HERE
        # #Do Bayesian block analysis
        # myblocks=bb.bblock(df,data_mode=3)
        # myblocks.find_blocks(fp_rate = fp_rate)
        # ########
        
        # #Convert rate from log of flux to flux
        # r_vec=10**myblocks.rate_vec
        # #Remove possible floating point erros by rounding r_vec down
        # r_vec=np.floor(r_vec)
        
        
        # c_points=myblocks.change_points
        
        
            
        
        
        
        # #Array for saving points for fitting a GPD
        # f_points=[]
        # #f_x=[]
        # #Pick samples only from blocks that are equal to or larger than u
        # for ii in range(len(c_points)-1):
        #     #Find maximum between the current change points and record it
        #     #if it is larger than the threshold u
        #     if r_vec[ii]>u:
        #         c_max=10**np.amax(df.FLUX[c_points[ii]:c_points[ii+1]])
        #         f_points.append(c_max)
                
                
        # f_points=np.asarray(f_points)
        
        #####UNCOMMENT ENDS HERE
        
        
        
        #Array for points from declustering
        fd_points=maxis
        # #f_x=[]
        # #Pick samples only from maxima that are equal to or larger than u
        # for ii in range(len(indis)-1):
        #     #Find maximum between the current declustering points and record it
        #     #if it is larger than the threshold u
        #     if maxis[ii]>u:
        #         c_max=10**np.amax(df.FLUX[indis[ii]:indis[ii+1]])
        #         fd_points.append(c_max)
                
                
        fd_points=np.asarray(fd_points)
        f_points=fd_points
        f_min_u=f_points-u
        
        
        #Convert log values to measurement values
        #and create plottting parameters
        stop=int(np.amax(f_min_u))
        start=int(np.amin(f_min_u))
        #start=int(10**np.amin(df.FLUX)-u)
              
        num=(stop-start)*10
        plot_x=np.linspace(start=start,stop=stop,num=num)
        plot_q=np.ones(len(df.TIME))*u
        
        #Plot the GPD data and fit a GPD to it
        #Calculate the autocorrelation
        if len(f_points)==0:
            autoc='nan'
            print("Length of f_points is 0.")
            #return
        else:
            autoc=sm.tsa.acf(f_min_u,nlags=1)
            autoc=autoc[1]
            #autoc=ce.my_autoc_lag(f_points,1)
            #autoc=autoc[0,1]
        
        # #Set x values for plotting
        # x_min=np.amin(plot_x)
        # x_max=np.amax(plot_x)
        # N=int((x_max-x_min)*10)
        # GPD_x=np.linspace(int(x_min),int(x_max),N)
        
        if len(f_points)>0:
            #Fit a GPD to the data
            shape,loc,scale=stats.genpareto.fit(f_min_u)
            #Asymptotic variance of the shape and scale parameters parameter
            var_s=(1+shape)**2/len(f_points)
            sd_s=np.sqrt(var_s)
            var_sc=(1+shape)*2*scale**2/len(f_points)
            sd_sc=np.sqrt(var_sc)
    
        #Set the correct xlbl for GPD and the ylabel (sic!) for the BB figures
        if flux_name=='fe130':
            xlbl='e flux (keV−1cm−1s−1ster−1)'
        elif flux_name=='fe1p2':
            xlbl='e flux (keV−1cm−1s−1ster−1)' 
        elif flux_name=='sgeo':
            xlbl='ULF spectral power for sgeo'
        elif flux_name=='sgr':
            xlbl='ULF spectral power for sgr'
    
        #while err in errs:
        if plt_GPD is not None:
            #Do plots of the GPD fit
            fig=plt.figure(11)
            ax=fig.add_subplot(1,1,1)
            #Plot data as a histogram
            ax.hist(f_min_u,bins='fd',density=True,label='experimental density')
            ax.set_xlim(start,stop)
            #Plot GPD base on the fit
            GPD_f=stats.genpareto.pdf(plot_x,shape,loc,scale)
            ax.set_xscale('log')
            ax.plot(plot_x,GPD_f,label='as per fit')
            #Add textbox
            n_uni=len(np.unique(f_points))
            textstr=f"shape=\
                    {np.round(shape,3)}+-{np.round(sd_s,5)}\n\
                    scale={np.round(scale)}+-{np.round(sd_sc,5)}\n\
                    N={np.round(len(f_points),0)}\n\
                    u={np.round(u,0)}\n\
                    unique={np.round(n_uni,0)}\n\
                    autoc_lag_1={np.round(autoc,3)}"
                    
                    #min steps={np.round(steps_min,0)}\n\
                    #max steps={np.round(steps_max,0)}\n\
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            # place a text box in upper left in axes coords
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                    verticalalignment='top', bbox=props)
                
            #Set figure layout
            ax.set_xlabel(xlbl)
            ax.set_ylabel('probability density')
            ax.set_xscale('log')
            ax.legend(loc='best')
            ax.set_xlim(np.amin(f_min_u),np.max(f_min_u))
            save_name=f"{flux_name}_"+t_add+f"_u={u}_err={err}_thr={thr}_{ns}sigma_{name_add}{avg_txt}{avg_d}_GPD_fit_fp={fp_rate}.pdf"
            plt.savefig(save_name)
            plt.close(11)
            
            #Create quantiles of the data and the fitted distribution
            
            
            #Plot a qq plot of the sample and the fitted GPD
            fig_n=plt.figure()
            ax_n=fig_n.add_subplot(1,1,1)
            #ss.probplot()f_min_u,
            # sm.qqplot(f_min_u,ssd.genpareto(shape),line='s',ax=ax) 
            # ax_n.set_xscale('log')
            # ax_n.set_xlabel('experimental quantiles of excesses')
            # ax_n.set_xlabel('theoretical quantiles of fitted GPD')
            # ax_n.set_title('qq plot, u={u}')
            
            save_name=f"qq_{flux_name}_"+t_add+f"_u={u}_err={err}_thr={thr}_{ns}sigma_{name_add}{avg_txt}{avg_d}_GPD_fit_fp={fp_rate}.pdf"
            fts.plot_my_qq_GPD(f_min_u,u,shape,scale,pname=save_name,logax=None)
            save_name=f"pp_{flux_name}_"+t_add+f"_u={u}_err={err}_thr={thr}_{ns}sigma_{name_add}{avg_txt}{avg_d}_GPD_fit_fp={fp_rate}.pdf"
            fts.plot_my_pp_GPD(f_min_u,shape,scale,pname=save_name,logax=None)
            # plt.savefig(save_name)
            
            # #Plot a probability plot of the sample and the fitted GPD
            # fig_n=plt.figure()
            # ax_n=fig_n.add_subplot(1,1,1)
            # sm.qqplot(f_min_u,ssd.genpareto,line='r',ax=ax_n,shape=shape,loc=loc,scale=scale) 
            # ax_n.set_xscale('log')
            # ax_n.set_xlabel('experimental quantiles of excesses')
            # ax_n.set_xlabel('theoretical quantiles of fitted GPD')
            # ax_n.set_title('qq plot, u={u}')
            # save_name=f"qq_{flux_name}_"+t_add+f"_u={u}_err={err}_thr={thr}_{ns}sigma_{name_add}{avg_txt}{avg_d}_GPD_fit_fp={fp_rate}.pdf"
            # plt.savefig(save_name)
            
    
        ####UNCOMMENT STARTS HERE
        #Create plotting vector from the rate vector r_vec
        # r_vec_plot=np.empty([0])
        # for ii in range(len(r_vec)):
        #     if ii==len(r_vec)-1:
        #         c_diff=1
        #     else:
        #         c_diff=c_points[ii+1]-c_points[ii]
        #     r_vec_plot=np.concatenate((r_vec_plot,r_vec[ii]*np.ones(c_diff)))
        ####UNCOMMENT ENDS HERE
            
        # #Do plotting vector for the declustered data
        # rd_vec_plot=np.empty([0])
        # for ii in range(len(maxis)):
        #     if ii==len(maxis)-1:
        #         c_diff=1
        #     else:
        #         c_diff=int(indis[ii+1])-int(indis[ii])
        #     rd_vec_plot=np.concatenate((rd_vec_plot,maxis[ii]*np.ones(c_diff)))
    
        
        curr_ax=ax12[kk,jj]
        curr_ax.scatter(df.TIME,10**df.FLUX,s=0.2,label='data')
        #curr_ax.plot(df.TIME,r_vec_plot,label='Bayesian blocks',color='r',linewidth=0.5)
        curr_ax.plot(df.TIME,plot_q,label='Threshold u',color='k',linewidth=0.5)
        
        xs=np.empty([0])
        for ii in range(len(indis)):
            #np.amax converts the esries type value into a ndarray compatible value
            add=np.amax(df.TIME[indis[ii]:indis[ii]+1])
            #add=int(add)
            xs=np.concatenate((xs,add*np.ones(1)))
        #xs=int(xs)
        #xs=df.loc[lambda df:df.iloc['TIME']==indis,:]
        curr_ax.scatter(xs,maxis,label=f"declustering ({thr} h)",color='orange',linewidth=0.5)
        
        
        #Set titles to subplots
        if flux_name=='fe130':     
            txt=f"iqr from {err} d"#", threshold value of u={np.round(u)}"
        elif flux_name=='fe1p2':
            txt=f"iqr from {err} d "#", threshold value of u={np.round(u)}"
        else:
            txt=f"iqr from {err} d"
        curr_ax.set_title(txt,fontsize=8)
        
        
        #Adjust subplots row and column
        jj=jj+1
        if jj==n_subpl:
            jj=0
            kk=kk+1
            
        #Collect garbage
        gc.collect()
         
    #######
    inds=np.linspace(0,n_subpl,n_subpl+1,dtype=int)
    med=int(np.floor(np.median(inds)))
    maxi=int(np.amax(inds)-1)
    ax12[med,0].set_ylabel(xlbl,labelpad=1.0) #sic!
    ax12[maxi,0].set_xlabel('time stamps (h from beginning)',labelpad=1.0)
    #Save the histograms
    fig12.suptitle(f"Data and Bayesian blocks {flux} for different errors", fontsize=12)
    #Add a common legend for the suplots
    handles, labels = ax12[med,0].get_legend_handles_labels()
    fig12.legend(handles, labels, loc='lower right')
    #Save the figure
    save_name=f"{flux_name}_"+t_add+f"_errs_thr={thr}_{ns}sigma_{name_add}{avg_txt}{avg_d}_BB_plots{cut_txt}_fp={fp_rate}.pdf"
    fig12.savefig(save_name)
    plt.close()
    #######
        
        #X values for plotting a horizontal line
        #print(myblocks.change_points)
        #Return number of change points, r_vec and data
        #N_c=len(c_points)
        #return fig12
    #Collect garbage
    gc.collect()

fp_rate=0.003
move_sd=1
mv_avg=None #If using mv_avg, set err to 42
slice_df=None
cut_first=None
plt_GPD=1
t_add=None
errs=np.linspace(20,20,1,dtype=int) #-2,0,9

flux="fe130"
path="/Users/mikkosavola/Documents/Opiskelu_HY/PhD/EVT_flux_project/"\
                    "data/pkl"
path.replace(" ", "")
path.replace("\n", "")
os.chdir(path)

#N_c=[]
thra=np.array([5,10,20])
for ii in range(len(thra)):
    thr=thra[ii]
    plot_GPD_bb(errs,flux,plt_GPD=plt_GPD,slice_df=slice_df,t_add=t_add,move_sd=move_sd,cut_first=cut_first,mv_avg=mv_avg,fp_rate=fp_rate,thr=thr)
    




    
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
#Fit a generalised Pareto distribution to the given data
#and get the distribution's parameters as a result
#q is a parameter form 0 to 1 to define what quantile is the lowest
#bound, and decl indicates whether the data is declustered or not

def plot_GPD_fit(q,decl):
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
    q=q
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
    fe130_orig=fe130
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
    
    if decl==1:
        #decluster fe130 data
        #set number of steps, within which the declustering is done
        decl_lim=4
        fe130_decl=np.empty(1)
        #limit_e=np.amax(fe130)+1 #For ensuring that find_min_dst works properly
        #Set the segment length for declustering
        #l_seg=len(fe130)
        #segments=int(np.ceil(l_seg/decl_lim))
        #Set n and k for find_Dst_min
        #n=24
        #k=11
        test=fe130[0]
        count=0
        for ii in range(len(fe130)):
            if fe130[ii]<test:
                count=count+1
                if count==decl_lim:
                    fe130_decl=np.append(fe130_decl,np.array(test))
                    if ii==len(fe130)-1:
                        test=fe130[ii]
                    else:
                        test=fe130[ii+1]
                    count=0
            elif fe130[ii]>=test:
                test=fe130[ii]
                count=0
        fe130_decl=fe130_decl[1:]
        fe130=fe130_decl
                    
# ############
#     ######FIND DST MAXIMA AND MINIMA AND THEIR INDICES IN THE ULF RECORD
        # for ii in range(segments-1):
        #     #Look for maximum e flux
        #     #Address the final segment separately, in case its length is less than n
        #     if ii==segments-2:
        #         ind_max=am.find_min_dst(-fe130[ii*n:l_seg],n,k,limit_e)
        #     else:
        #         ind_max=am.find_min_dst(-fe130[ii*n:(ii+2)*n],n,k,limit_e)
        #     if ind_max is False:
        #         continue
        #     else:
        #         ind_max=ind_max+ii*n
        #         fe130_decl=np.append(fe130_decl,np.array(fe130[ind_max]))
# ##############               
        #fe130_decl=np.asarray(fe130_decl)
        
    
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
    # dual_sgeo=fts.dual_dist(sgeo)
    # dual_sgr=fts.dual_dist(sgr)
    # dual_fe1p2=fts.dual_dist(fe1p2)
    # dual_fe130=fts.dual_dist(fe130)
    
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
    

    
    ##Enter here the values to be plotted and fit as a GPD
    
    #xi and sd are from a bootstrap test
    xi=0.297
    sd=0.049
    plot_x=fe130 #fe130, fe1p2, sgeo or sgr
    plot_x_orig=fe130_orig
    smpl='fe130'
    #linewidth
    lw=0.5
    #Number of standard deviations
    n_sd=1
    #labels for plot
    xlbl='e flux (keV−1cm−1s−1ster−1)'
    
    #Find indices of the remaining data points in the original array
    #to see how far apart are the closest point
    steps_min='nan'
    steps_max='nan'
    if decl==0:
        ind_test=np.zeros(len(plot_x_orig))
        ind_test[plot_x_orig>=np.amin(plot_x)]=1
        steps_min=1e5
        steps_max=0
        #Loop through all test values
        for ii in range(len(ind_test)):
            if ind_test[ii]==1:
               #Check when the next '1' is found in the sequence
                for kk in range(ii+1,len(ind_test)):
                    if ind_test[kk]==1:
                        steps_new=kk-ii
                        #Evaluate the new smallest and largest gap between '1's
                        if steps_new<steps_min:
                            steps_min=steps_new
                        if steps_new>steps_max:
                            steps_max=steps_new
    
    
    #Calculate the autocorrelation and pick its maximum absolute value
    autoc=sm.tsa.acf(plot_x,nlags=1)
    
    
    #Set x values for plotting
    x_min=np.amin(plot_x)
    x_max=np.amax(plot_x)
    N=int((x_max-x_min)*10)
    GPD_x=np.linspace(int(x_min),int(x_max),N)
    
    #Fit a GPD tp the data
    shape,loc,scale=stats.genpareto.fit(plot_x)
    #Asymptotic variance of the shape and scale parameters parameter
    var_s=(1+shape)**2/len(plot_x)
    sd_s=np.sqrt(var_s)
    var_sc=(1+shape)*2/len(plot_x)
    
    
    #Do plots
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    #Plot data as a histogram
    ax.hist(plot_x,bins='fd',density=True,label='experimental density')
    #Plot GPD as per bootstrap results
    GPD_d=stats.genpareto.pdf(plot_x,xi,loc=x_min)
    ax.set_xlim(x_min,10**7)
    #ax.plot(plot_x,GPD_d,label='as per bootstrap',linewidth=lw)
    #ax.plot(plot_x,GPD_d-n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    #ax.plot(plot_x,GPD_d+n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    #Plot GPD base on the fit
    GPD_f=stats.genpareto.pdf(GPD_x,shape,loc,scale)
    ax.set_xscale('log')
    ax.plot(GPD_x,GPD_f,label='as per fit')
    #Add textbox
    textstr=f"shape=\
            {np.round(shape,3)}+-{np.round(sd_s,5)}\n\
            scale={np.round(scale)}+-{np.round(sd_s,5)}\n\
            q={np.round(q,6)}\n\
            N={np.round(len(plot_x),0)}\n\
            u={np.round(np.amin(plot_x,0))}\n\
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
    ax.legend()
    if decl==1:
        plt.savefig(f"{smpl}_GPD_fit_q={np.round(q,6)}_decl_lim={decl_lim}.pdf")
    if decl==0:
        plt.savefig(f"{smpl}_GPD_fit_q={np.round(q,6)}.pdf")
    plt.close()
    
decl=1
# q=0.90
# while q<0.99:
#     plot_GPD_fit(q,decl)
#     q=q+0.1
q=0.98
while q<0.999:
    plot_GPD_fit(q,decl)
    q=q+0.01
q=0.995
while q<0.9999:
    plot_GPD_fit(q,decl)
    q=q+0.001
q=0.999
while q<0.9999:
    plot_GPD_fit(q,decl)
    q=q+0.0001
    
    




    
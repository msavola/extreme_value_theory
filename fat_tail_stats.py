#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:41:23 2021

@author: mikkosavola
"""

import aux_mi as am
import main_mi as mm
import matplotlib.pyplot as plt
#import main_mi as mm
import numpy as np
import numpy.random as nr
import scipy.stats as st
import scipy.stats.distributions as ssd
import os
import statsmodels.api as sm
from statsmodels.base.model import GenericLikelihoodModel
import scipy.special as ss
import scipy.optimize as sop
#This module contains functions for approximating whether a
#distribution underlying a time series is fat-tailed

import time
import functools
from typing import Callable

# Code sample  for 'timer' is from here
# https://realpython.com/primer-on-python-decorators/
#This function times the execution of another function
def timer(func: Callable) -> Callable:
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start = time.perf_counter()
        print(f"Execution of {func.__name__!r} started at {start}")
        
        result = func(*args, **kwargs)
        
        end = time.perf_counter()
        print(f"\tEnded at {end}")
        run_time = end-start
        print(f"\tTime taken {run_time:.4f} seconds.")
        return result
    return wrapper_timer

# Here another decorator for counting calculation iteration?


@timer
def clean_data(x,a,big):
    #Removes nans and too big values and returns the data
    #and the year stamps
    
    #use a dummy object to access object functions
    dmm=mm.My_Mi('1',1,1,'2',2,2)
    #Remove nans and too big values
    x=dmm._big2nans(x,big)
    arr=dmm._remnans_n(x,a)
    x=arr[:,0]
    a=arr[:,1]
    
    #Return cleaned data and the year labels
    return x,a

def add_noise(x,lim):
    #Adds small random noise of size lim to x
    
    #Seed the PRNG
    sq1=nr.SeedSequence()
    #sq2=nr.seedsequence(sq1.entropy)
    seed=sq1.generate_state(1)
    rng=nr.default_rng(seed=seed)
    lx=len(x)
    #Create small noise
    noise=rng.uniform(low=0,high=1,size=lx)
    for ii in range(lx):
        noise[ii]=scale_noise(x[ii],noise[ii],lim)
    x=x+noise
    
    return x


def scale_noise(x,n,lim):
    #Scales the noise added by function add_noise to be
    #proportional to the original data
    scale=n/x
    scale=adj_scale(scale,lim)
    return scale
  
  
def adj_scale(scale,lim):
    #Check whether scale is between lim and lim/10 and
    #adjust accordingly
    status=False
    u_limit=lim
    l_limit=u_limit/10
    while status==False:
        if abs(scale)>u_limit:
            scale=scale/10
        elif abs(scale)<l_limit:
            scale=scale*10
        else:
            status=True
    return scale
    

@timer
def dual_dist(x,H=None):
    #Produces the unbounded a dual distribution of x
    #It's assumed that L is tha lower bound of X and
    #H is the upper finite bound of x
    L=min(x)
    if H==None:
        H=max(x)
    y=L-H*np.log((H-x)/(H-L))
    return y

def mef(x,u):
    #Mean excess function of values of x
    #u is the threshold value
    x_pos=x
    x_pos=x_pos[x_pos>u]
    den=len(x_pos)
    num=np.sum(x_pos)-den*u  
    if num==0 and den==0:
        e_u=0
    else:
        e_u=num/den
        
    #Convert array to scalar    
    if type(e_u)==np.ndarray:
        e_u=e_u[0]
    return e_u

def m2s(x,p):
    #Maximum to sum
    mp=np.amax(x**p)
    sp=np.sum(x**p)
    #Handle division by zero
    if mp==0 and sp==0:
        rp=0
    elif sp==0:
        rp=float('inf')
        print('########\nInifinity encountered in m2s plot\n#######')
    else:
        rp=mp/sp
    
    return rp

def zipf(x,l):
    #Returns the empirical probability of the survival function
    scale=len(x)
    #prob=sum(np.where(x>l,1,0))/scale
    prob=np.zeros([scale,1])
    prob[x>l]=1
    prob=np.sum(prob)/scale
    return prob
   

# def qq_exp(x):
#     #Calculates the empirical quantiles of x and the corresponding
#     #exponential quantiles
#     lx=len(x)
#     #Avoid infinity in exp by adding/subtracting eps
#     #from quantile limits 
#     eps=10**(-7)
#     qs=np.linspace(0+eps,1-eps,lx)
#     quants=np.quantile(x,qs,interpolation='linear')
        
#     #Compute the theoretical quantiles
#     q_exp=-np.log(1-qs)  
    
#     return quants,q_exp

@timer
def qq_gauss_plot(x,fname):
    #Calculates the empirical quantiles of x and the corresponding
    #normalised Gaussian quantiles
    
    #Flatten x to make it a 1D array
    x=np.ndarray.flatten(x)
    
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    #Fix the location and scale of the theoretical dist
    #ribution and also add a scaled "45 degree" line
    sm.qqplot(x,ssd.norm,line='s',ax=ax,loc=0,scale=1)
    
    #plt.xscale('log')
    ax.set_title('qq plot, w/ theoretical quantiles from N(0,1)')
    fig.savefig('qq_gauss_'+fname+'.pdf')
    plt.close(fig)
    #fig.show()
    
@timer
def qq_exp_plot(x,fname):
    #Calculates the empirical quantiles and the corresponding
    #exponential quantiles
    
    #Flatten x to make it a 1D array
    x=np.ndarray.flatten(x)
    #Create figure for plotting
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    #Fix the location and scale of the theoretical dist
    #ribution and also add a scaled "45 degree" line
    sm.qqplot(x,ssd.expon,line='s',ax=ax,loc=0,scale=1)
    
    #plt.xscale('log')
    ax.set_title('qq plot, w/ theoretical quantiles from exp(-x)')
    fig.savefig('qq_exp_'+fname+'.pdf')
    plt.close(fig)
    #fig.show()
    
def my_pp_GPD(exc,shape,scale):
    #Create probabilities of a sample and the fitted GPD
    #exc are the exceedances above the threshold u
    #The exceedances are sorted into an increasing order
    smp=np.sort(exc)
    p_smp=np.empty([0])
    p_fit=np.empty([0])
    k=len(exc)
    
    for ii in range(k):
        y=smp[ii]
        p=(ii+1)/(k+1)
        Gi=1-(1+shape*y/scale)**(-1/shape)
        p_smp=np.concatenate((p_smp,np.ones(1)*p))
        p_fit=np.concatenate((p_fit,np.ones(1)*Gi))
    
    return p_fit,p_smp
    
def my_qq_GPD(exc,u,shape,scale):
    #Create quantiles of a sample and the fitted GPD
    #exc are the exceedances above the threshold u
    #The exceedances are sorted into an increasing order
    
    q_smp=np.sort(exc)+u
    q_fit=np.empty([0])
    k=len(exc)
    
    #Parameters for quantile estimation
    #a=b=1/3 yeilds approximately median-unbiased estimates
    #for any distribution see Hyndman&Fan 1996
    a=1/3
    b=1/3
    
    for ii in range(k):
        y=(ii+1-a)/(k+1-a-b)
        Gi=u+scale/shape*((1-y)**(-shape)-1)
        q_fit=np.concatenate((q_fit,np.ones(1)*Gi))
    
    return q_fit,q_smp

def my_regress(x,y,log=None):
    if log==None:
        #Do a linear regression
        slp,icp,r,p,se=st.linregress(x,y)
        #Usable for Scipy >1.7, instead of the above lines
        # res=st.linregress(q_fit,q_smp,alternative='greater')
        # slp=res.slope
        # icp=res.intercept
        # r=res.rvalue
        # p=res.value
        # slp_se=res.stderr
        # icp_se=res.intercept_stderr
    if log=='x':
        slp,icp,r,p,se=st.linregress(np.log10(x),y)
    elif log=='y':
        slp,icp,r,p,se=st.linregress(x,np.log10(y))
    elif log=='yx':
        slp,icp,r,p,se=st.linregress(np.log10(x),np.log10(y))
    elif log=='xy':
        slp,icp,r,p,se=st.linregress(np.log10(x),np.log10(y))
    
    return slp,icp,r,p,se

def plot_my_qq_GPD(exc,u,shape,scale,pname=None,logax=None):
    #Plot the quantiles and fit a regression line
    q_fit,q_smp=my_qq_GPD(exc,u,shape,scale)
    
    #Do a linear regression
    slp,icp,r,p,se=my_regress(q_fit,q_smp,logax)
    R2=r**2
    
    #Values for plotting the regression line
    start=int(np.amin(q_fit))
    stop=int(np.amax(q_fit))
    num=int((stop-start)*10)
    x_plot=np.linspace(start,stop,num)
    y_plot=slp*x_plot+icp
    
    ###This is still to be refined, the case logax=='x' seems to work though
    # if logax=='x':
    #     y_plot=slp*np.log10(x_plot)+icp
    # if logax=='y':
    #     y_plot=np.log10(y_plot)
    
    
    
    #Plot the results
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    ax.scatter(q_fit,q_smp,color='black')
    ax.plot(x_plot,y_plot,color='blue')
    ax.set_title('Quantile-quantile plot (median-unbiased)')
    ax.set_xlabel('model quantiles')
    ax.set_ylabel('sample quantiles')
    
    #Adjust axis, if required, and do the regression accordingly
    if logax=='x':
        ax.set_xscale('log')
    elif logax=='y':
        ax.set_yscale('log')
    elif logax=='xy':
        ax.set_xscale('log')
        ax.set_yscale('log')
    elif logax=='yx':
        ax.set_xscale('log')
        ax.set_yscale('log')
                 
        
    #Add legend   
    plt.legend(loc='best')
    
    #Add textbox
    textstr=f"R^2={np.round(R2,3)}\
              \nslope of regr.={np.round(slp,3)}  "
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    
    if pname==None:
        pname='newfig.pdf'
    
    plt.savefig(pname)
    return

def plot_my_pp_GPD(exc,shape,scale,pname=None,logax=None):
    #Plot the probabilities and fit a regression line
    p_fit,p_smp=my_pp_GPD(exc,shape,scale)
    
    #Do a linear regression
    slp,icp,r,p,se=my_regress(p_fit,p_smp,logax)
    R2=r**2
    
    #Values for plotting the regression line
    x_plot=np.linspace(0,1,100)
    y_plot=slp*x_plot+icp
    
    ###This is still to be refined, the case logax=='x' seems to work though
    # if logax=='x':
    #     y_plot=slp*np.log10(x_plot)+icp
    # if logax=='y':
    #     y_plot=np.log10(y_plot)
    
    
    
    #Plot the results
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    ax.scatter(p_fit,p_smp,color='black')
    ax.plot(x_plot,y_plot,color='blue')
    ax.set_title('Probability plot')
    ax.set_xlabel('model probabilities')
    ax.set_ylabel('sample probabilities')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    
    #Adjust axis, if required, and do the regression accordingly
    if logax=='x':
        ax.set_xscale('log')
    elif logax=='y':
        ax.set_yscale('log')
    elif logax=='xy':
        ax.set_xscale('log')
        ax.set_yscale('log')
    elif logax=='yx':
        ax.set_xscale('log')
        ax.set_yscale('log')
                 
        
    #Add legend   
    plt.legend(loc='best')
    
    #Add textbox
    textstr=f"R^2={np.round(R2,3)}\
              \nslope of regr.={np.round(slp,3)}  "
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    
    if pname==None:
        pname='newfig.pdf'
    
    plt.savefig(pname)
    return
    
def find_large_axis(x):
    #Finds and returns the axis that has the largesrt size
    s=np.shape(x)
    res=np.argmax(s)
    return res
    
def asc_sort(x):
    #Sort the elements of array x in ascending order according
    #to the axis with the largestt size
    ax=find_large_axis(x)
    res=np.sort(x,axis=ax)
    return res
    
def desc_sort(x):
    #Sort the elements of array x in descending order according
    #to the axis with the largest size
    ax=find_large_axis(x)
    res=np.flip(np.sort(x,axis=ax))
    return res

def hill_est(x,k):
    #Returns the Hill estimator for the
    #parameter xi of a Generalised Pareto distribution
    
    #Check that 0<=k<=n-1
    if k<0:
        print('k must be between 1 and sample size')
        return
    if k>len(x)-1:
        print('k must be between 1 and sample size')
        return x
    
    x_d=desc_sort(x)
    xi=0
    #Add 1 to range of k to have the last term in the sum included
    for ii in range(k+1):
        xi=xi+np.log(x_d[ii])-np.log(x_d[k])
    xi=xi/k
    return xi

@timer
def fd_hist(x,ttl,units):
    #Plot a histogram of x using the Freedman-Diaconis rule for binning
    fig=plt.figure()
    
    ax=fig.add_subplot(1, 1, 1)
    plt.hist(x,bins='fd',density=True)
    #plt.hist(x,density=True)
    ax.set_title('Normalised histogram (max value at end of x axis)')
    x_txt='bin values ('+units+')'
    y_txt='normalised frequency'
    ax.set_xlabel(x_txt)
    ax.set_ylabel(y_txt)
    x_max=np.amax(x)
    x_min=np.amin(x)
    ax.set_xlim(x_min,x_max)
    #Set log scale, if required
    # if x_scale=='log':
    #     ax.set_xscale('log')
    #     ttl=ttl+'_xlog_'
    # if y_scale=='log':
    #     tt=ttl+'_ylog_'
    #     ax.set_yscale('log') 
    fig.savefig('hist_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    return
    
def pick_est(x,k):
    #Returns the Pickand's estimator for the
    #parameter xi of a Generalised Pareto distribution
    
    #Re-index k to correctly apply the formula for Pickand's estimator
    #k=k+1
    
    #Check that k has a valid value
    if k<1:
        print('k must be between 1 and 1/4 of sample size')
        return
    if k>len(x)/4:
        print('k must be between 1 and 1/4 of sample size')
        import pdb
        pdb.set_trace()
        return
    
    #Evaluate the estimator's value
    x_a=desc_sort(x)
    num=x_a[k-1]-x_a[2*k-1]
    den=x_a[2*k-1]-x_a[4*k-1]
    xi=1/np.log(2)*np.log(num/den)
    return xi

def dedh_est(x,k):
    #Returns the Dekkers'Einmahl-de Haan estimator for the
    #parameter xi of a Generalised Pareto distribution
    
    
    #Evaluate the estimator's value
    x_a=desc_sort(x)
    
    H1=0
    H2=0
    for ii in range(k):
        H1=H1+1/k*(np.log(x_a[ii])-np.log(x_a[k+1]))
    for jj in range(k):
        H2=H2+1/k*(np.log(x_a[jj])-np.log(x_a[k+1]))**2
    xi=1+H1+0.5*(H1**2/H2-1)**(-1)
    return xi

# def qq_exp_plot(x,ttl):
#     lx=len(x)
#     #Calculate the quantiles
#     #qx=np.empty(lx)
#     #qy=np.empty(lx)
#     qx,qy=qq_exp(x)
#     #Create plotting points for a 45 degree line
#     yp_max=np.amax(x)
#     yp_min=np.amin(x)
#     yp=np.linspace(yp_min,yp_max,lx)
#     #Plot the quantiles with a 45 degree line
#     fig=plt.figure()
#     ax=fig.add_subplot(1, 1, 1)
#     ax.scatter(qy,qx)
#     ax.plot(yp,yp,color='red',label='45 degree line')
#     ax.set_title('qq plot')
#     ax.set_xlabel('theoretical quantiles of exp(-x)')
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     ax.set_ylabel('empirical quantiles')
#     fig.savefig('qq_exp_'+ttl+'.pdf')
#     fig.show()
#     return

# def qq_gauss_plot(x,ttl):
#     lx=len(x)
#     #Calculate the quantiles
#     #qx=np.empty(lx)
#     #qy=np.empty(lx)
#     qx,qy=qq_gauss(x)
#     #Create plotting points for a 45 degree line
#     yp_max=np.amax(x)
#     yp_min=np.amin(x)
#     yp=np.linspace(yp_min,yp_max,lx)
#     #Plot the quantiles with a 45 degree line
#     fig=plt.figure()
#     ax=fig.add_subplot(1, 1, 1)
#     ax.scatter(qy,qx)
#     ax.plot(yp,yp,color='red',label='45 degree line')
#     ax.set_title('qq plot')
#     ax.set_xlabel('theoretical quantiles of N(0,1)')
#     #ax.set_xscale('log')
#     #ax.set_yscale('log')
#     ax.set_ylabel('empirical quantiles')
#     fig.savefig('qq_gauss_'+ttl+'.pdf')
#     fig.show()
#     return

@timer
def me_plot(x,ttl,x_lbl):
    xis=[]
    x_sort=asc_sort(x)
    for ii in range(len(x)):
        xis.append(mef(x_sort,x_sort[ii]))  
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.scatter(x_sort,xis)
    ax.set_title('me-plot')
    x_txt='threshold ('+x_lbl+')'
    y_txt='mean excess'
    ax.set_xlabel(x_txt)
    ax.set_ylabel(y_txt)
    #Set log scale, if required
    # if x_scale=='log':
    #     ax.set_xscale('log')
    #     ttl=ttl+'_xlog_'
    # if y_scale=='log':
    #     tt=ttl+'_ylog_'
    #     ax.set_yscale('log') 
    fig.savefig('me_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    return

@timer
def me_plot_ll(x,ttl,x_lbl):
    xis=[]
    x_sort=asc_sort(x)
    for ii in range(len(x)):
        xis.append(mef(x_sort,x_sort[ii]))
    #xis=np.array([xis])  
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.scatter(x_sort,xis)
    ax.set_title('me-plot')
    x_txt='threshold in log scale ('+x_lbl+')'
    y_txt='mean excess'
    ax.set_xlabel(x_txt)
    ax.set_ylabel(y_txt)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    fig.savefig('me_log_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    return

@timer
def m2s_plot(x,p,ttl,n_x=None):
    m2s_s=[]
    x=np.array([x])
    #x=asc_sort(x)
    
    #Different options for n_x
    if n_x is None:
        n_x=max(np.shape(x))
    else:
        n_x=n_x
    
    x_p=np.linspace(1,n_x,n_x)
    for ii in range(n_x):
        m2s_s.append(m2s(x[0,0:ii+1],p))
    m2s_s=np.array([m2s_s])    
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,m2s_s[0,:])
    ax.set_title('ms-plot for p='+str(p))
    ax.set_xlabel('n')
    ax.set_ylabel('R_n')
    fig.savefig('m2s_'+ttl+'p='+str(p)+'.pdf')
    plt.close(fig)
    #fig.show()
    return
   
@timer 
def zipf_plot(x,xlbl,ttl):
    zipfs=[]
    for ii in range(len(x)):
        zipfs.append(zipf(x,x[ii]))
        #print('zipf',ii)
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.scatter(x,zipfs)
    ax.set_title('zipf-plot')
    ax.set_xlabel(xlbl)
    ax.set_ylabel('empirical survival rate') 
    fig.savefig('zipf_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
  
@timer
def zipf_plot_ll(x,xlbl,ttl):
    zipfs=[]
    for ii in range(len(x)):
        zipfs.append(zipf(x,x[ii]))
        #print('zipf',ii)
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.scatter(x,zipfs)
    ax.set_title('zipf-plot')
    ax.set_xlabel(xlbl+', log scale')
    ax.set_ylabel('empirical survival rate, log scale') 
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig('zipf_log_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
 
@timer
def pick_plot(x,k,ttl):   
#Plots Pickand's estimator for xi  
    picks=[]
    var=[]
    for ii in range(1,k+1):
        picks.append(pick_est(x,ii))
        xi=picks[ii-1]
        num=xi**2*(2**(2*xi+1)+1)
        den=(2*(2**xi-1)*np.log(2))**2
        vari=num/den/np.sqrt(ii)
        var.append(vari)
    x_p=np.linspace(1,k,k)    
    fig=plt.figure()
    n_sd=2
    #Convert to arrays for calculations
    picks=np.asarray(picks)
    var=np.asarray(var)
    sd=var**0.5
    
    lw=0.5
    ylims=2.0
    
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,picks,label='estimate',linewidth=lw)
    ax.plot(x_p,picks-n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    ax.plot(x_p,picks+n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    ax.set_title('Plot of Pickand''s estimator for xi')
    ax.set_xlabel('k')
    ax.set_ylabel('xi')
    ax.set_ylim([-ylims,ylims])
    ax.legend()
    fig.savefig('pick_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()

@timer
def dedh_plot(x,k,ttl):   
#Plots Dekkers-Einhal-de Haan estimator for xi  
    dedhs=[]
    var=[]
    for ii in range(1,k+1):
        dedhs.append(dedh_est(x,ii))
        xi=dedhs[ii-1]
        if xi>=0:
            vari=(1+xi**2)/np.sqrt(ii)
        else:
            b=4-8*(1-2*xi)/(1-3*xi)+(5-11*xi)*(1-2*xi)/((1-3*xi)*(1-4*xi))
            vari=(1-xi)**2*(1-2*xi)*b/np.sqrt(ii)
        var.append(vari)
    x_p=np.linspace(1,k,k)    
    fig=plt.figure()
    n_sd=2
    #Convert to arrays for calculations
    dedhs=np.asarray(dedhs)
    var=np.asarray(var)
    sd=var**0.5
    
    lw=0.5
    ylims=2.0
    
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,dedhs,label='estimate',linewidth=lw)
    ax.plot(x_p,dedhs-n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    ax.plot(x_p,dedhs+n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    ax.set_title('Plot of DEdH''s estimator for xi')
    ax.set_xlabel('k')
    ax.set_ylabel('xi')
    ax.set_ylim([-ylims,ylims])
    ax.legend()
    fig.savefig('dedh_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    
    
@timer
def pick_plot_N(dist,a,scale,size,k,ttl,N,conf):   
#Plots Pickand's estimator for xi of a GPD (given at input dist)
#as the average of N runs and the conf % (converted to percentile limits high_conf
#and low_conf experimental confidence interval
    picks=np.empty([k,N])
    for jj in range(N):
        x=dist(a,scale=scale,size=size)
        for ii in range(1,k+1):
            picks[ii-1,jj]=pick_est(x,ii)
    x_p=np.linspace(1,k,k)    
    #Calculate the average values of xi and the conf % confidence interval
    avg_xi=np.mean(picks,axis=1)
    #Do a descending ordering of the xi values over N for all k
    picks_asc=picks
    for ii in range(k):
        picks_asc[ii,:]=asc_sort(picks[ii,:])
    #Plot the figure
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,avg_xi,label='avg. of xi')
    #Plot confidence intervals
    low_conf=int(np.ceil(((1-conf)/2)*N-1))
    high_conf=int(np.floor((conf/2+0.5)*N))
    low_per=int(low_conf/(N-1)*100)
    high_per=int(high_conf/(N-1)*100)
    #true_conf=(high_conf-low_conf)*100
    ax.plot(picks_asc[:,low_conf],label=str(high_per)+' percentile')
    ax.plot(picks_asc[:,high_conf],label=str(low_per)+' percentile')
    
    ax.set_title('Plot of average Pickand''s estimator for xi, N='+str(len(x))+', runs='+str(N))
    ax.set_xlabel('k')
    ax.set_ylabel('xi')
    ax.legend()
    fig.savefig('pick_N='+str(len(x))+'_runs='+str(N)+'_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
      
   
@timer   
def hill_plot(x,k,ttl):  
#Plots Hill's estimator for xi
    hills=[]
    var=[]
    for ii in range(k-1):
        hills.append(hill_est(x,ii+2)) 
        var.append(hills[ii]**2/np.sqrt(ii+1))
    x_p=np.linspace(2,k,k-1)   
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    n_sd=2
    #Convert to arrays for calculations
    hills=np.asarray(hills)
    var=np.asarray(var)
    sd=var**0.5
    
    lw=0.5
    
    ax.plot(x_p,hills,label='estimate',linewidth=lw)
    ax.plot(x_p,hills-n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    ax.plot(x_p,hills+n_sd*sd,linestyle='--',linewidth=lw,label=f"{n_sd} sigma CI")
    ax.set_title('Plot of Hill''s estimator for xi')
    ax.set_xlabel('k')
    ax.set_ylabel('xi')
    ax.legend()
    fig.savefig('hill_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    
@timer   
def alt_hill_plot(x,k,ttl):  
#Plots Hill's estimator for xi with a logarithmic axis for k
    hills=[]
    #th=np.empty([k,1])
    for ii in range(k-1):
        hills.append(hill_est(x,ii+2)) 
        #th[ii,0]=ii
    x_p=np.linspace(2,k,k-1)
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.semilogx(x_p,hills,base=10)
    ax.set_title('AltHill plot of Hill''s estimator for xi')
    ax.set_xlabel('k (log10 scale)')
    ax.set_ylabel('xi')
    fig.savefig('alt_hill_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    
@timer   
def alt_pick_plot(x,k,ttl):  
#Plots Hill's estimator for xi with a logarithmic axis for k
    picks=[]
    for ii in range(k):
        picks.append(pick_est(x,ii)) 
    x_p=np.linspace(1,k,k)
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.semilogx(x_p,picks,base=10)
    ax.set_title('AltPickand plot of Pickand''s estimator for xi')
    ax.set_xlabel('k (log10 scale)')
    ax.set_ylabel('xi')
    fig.savefig('alt_pick_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    
@timer   
def hill_plot_N(dist,a,scale,size,k,ttl,N,conf):   
#Plots Hills's estimator for xi of a GPD (given at input dist)
#as the average of N runs and the conf % (converted to percentile limits high_conf
#and low_conf experimental confidence interval
    hills=np.empty([k-1,N])
    for jj in range(N):
        x=dist(a=a,scale=scale,size=size)
        for ii in range(k-1):
            hills[ii,jj]=hill_est(x,ii+2)
    x_p=np.linspace(2,k,k-1)    
    #Calculate the average values of xi and the conf % confidence interval
    avg_xi=np.mean(hills,axis=1)
    #Do a descending ordering of the xi values over N for all k
    hills_asc=hills
    for ii in range(k-1):
        hills_asc[ii,:]=asc_sort(hills[ii,:])
    #Plot the figure
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,avg_xi,label='avg. of xi')
    #Plot confidence intervals
    low_conf=int(np.ceil(((1-conf)/2)*N-1))
    high_conf=int(np.floor((conf/2+0.5)*N))
    low_per=int(low_conf/(N-1)*100)
    high_per=int(high_conf/(N-1)*100)
    #true_conf=(high_conf-low_conf)*100
    ax.plot(hills_asc[:,low_conf],label=str(high_per)+' percentile')
    ax.plot(hills_asc[:,high_conf],label=str(low_per)+' percentile')
    
    ax.set_title('Plot of average Hill''s estimator for xi, N='+str(len(x))+', runs='+str(N))
    ax.set_xlabel('k')
    ax.set_ylabel('xi')
    ax.legend()
    fig.savefig('hill_N='+str(len(x))+'_runs='+str(N)+'_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()

#Solve xi using MLE described in the book
#Modelling Extreme Events for Insurance and Finance (1997)
#p. 356
#y contains the excesses being inserted as observational data
def xi_MLE_est(x,k,y):
    #return xi and tau
    
    return [x[0]-1/k*np.sum(np.log(1-x[1]*y)),
            1/x[1]+1/k*(1/x[0]+1)*np.sum(y/(1-x[1]*y))]

@timer
#Plot xi as solved using MLE, with an initial guess:
#"guess" is the guess and y contains the excesses. k is the summing index
def xi_MLE_plot(guess,k,y,ttl):
    mles=np.empty([k,2])
    var=np.empty([k,2])
    for ii in range(k):
        xi,tau=sop.fsolve(xi_MLE_est,guess,args=(ii+1,y[0:ii+1])) 
        beta=-tau/xi
        mles[ii,:]=np.array([xi,beta])
        var[ii,:]=1/np.sqrt(ii+1)*np.array([(1+mles[ii,0])**2,1+mles[ii,0]*2])
    sd=var**0.5
    x_p=np.linspace(1,k,k)
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,mles[:,0],label='xi estimate values')
    n_sd=2
    ax.plot(x_p,mles[:,0]-n_sd*sd[:,0],color='blue',linestyle='--',label=f"{n_sd} s.d. CI for xi")
    ax.plot(x_p,mles[:,0]+n_sd*sd[:,0],color='blue',linestyle='--',label=f"{n_sd} s.d. CI for xi")
    ax.set_title('MLE estimate plot for xi')
    ax.set_xlabel('k')
    ax.set_ylabel('xi')
    #ylim=1.0
    #ax.set_ylim(-ylim,ylim)
    ax.legend()
    print(f"beta={np.round(mles[:,1],3)}")
    fig.savefig('mle_xi_'+ttl+'.pdf')
    plt.close(fig)
    
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(x_p,mles[:,1],label='sigma estimate values')
    ax.plot(x_p,mles[:,1]-n_sd*sd[:,1],color='red',linestyle='-.',label=f"{n_sd} s.d. CI for sigma")
    ax.plot(x_p,mles[:,1]+n_sd*sd[:,1],color='red',linestyle='-.',label=f"{n_sd} s.d. CI for sigma")
    ax.set_title('MLE estimate plot for sigma')
    ax.set_xlabel('k')
    ax.set_ylabel('sigma')
    ax.legend()
    fig.savefig('mle_sigma_'+ttl+'.pdf')
    plt.close(fig)
    #fig.show()
    
def do_plots(x,k,ttl,zipf_xlbl,m2s_p,a=None,p_year=None,m2s_n_x=None):  
    print('\n')
    print('Doing plots for',ttl)
    #Find indices of the entries corresponding to the  chosen years
    ind=np.empty([1,1],dtype=np.int32)
    if p_year is not None:
        if type(p_year)==np.ndarray:
            for ii in range(len(p_year)):
                ind=np.concatenate((ind,np.argwhere(a==p_year[ii])),axis=0)
            #Remove the empty array entry from the beginning
            ind=ind[1:,:]
        else:
            ind=np.argwhere(a==p_year)
    if a is not None:
        x_ind=x[ind]
        print('In year',str(p_year),'there are',str(len(x_ind)),'data points')
    else:
        x_ind=x   
        print('In the sample there are',str(len(x_ind)),'data points')
    print('\n')
    #Add sample size to ttl
    ttl=ttl+'_N='+str(len(x_ind))
    #Do histogram
    fd_hist(x_ind,ttl,zipf_xlbl)
    #Do me-plots
    me_plot(x_ind,ttl,zipf_xlbl)
    me_plot_ll(x_ind,ttl,zipf_xlbl)
    #Do quantile-quantile plots
    qq_gauss_plot(x_ind,ttl)
    qq_exp_plot(x_ind,ttl)
    #Do Zipf plots
    zipf_plot(x_ind,zipf_xlbl,ttl)
    zipf_plot_ll(x_ind,zipf_xlbl,ttl)
    #Do maximum to sum plots
    for ii in m2s_p:
        m2s_plot(x_ind,ii,ttl,m2s_n_x)
    #Do Pickand's plot and Hill plot
    pick_plot(x_ind,k,ttl)
    hill_plot(x_ind,k,ttl)
    #Do AltHill plot and AltPick plot
    alt_hill_plot(x_ind,k,ttl)
    alt_pick_plot(x_ind,k,ttl)
    #Do DEdH plot
    dedh_plot(x_ind,k,ttl)
    

def fit_GPD(x):
    #Fit a generalised Pareto distribution to the given data
    #and get the distribution's parameters as a result
    shape=0
    loc=0
    scale=0
        
    shape,loc,scale=ss.genpareto.fit(x)
    
    return scale,loc,scale
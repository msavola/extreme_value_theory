#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 10:59:15 2022

@author: mikkosavola
"""

import numpy as np

#Autocorrelation (normalised) of a time series x with itself for lag
#and for processes that are weak-sense stationary
def my_autoc_lag(x,lag):
    lx=len(x)
    
    x1=x[:lx-lag]
    x2=x[lag:]
    
    autoc=np.corrcoef(x1,x2)
    
    return autoc

#Do declustering of the data to pick out independent maxima above a threshold
#lim. Pick values that are local maxima cnt steps onward (and by consturcion
#also backward) and end a strictly increasing sequence.
#Pick the last potential maximum as a
#maximum when the array ends.
def declustering(x,lim,cnt):
    lx=np.amax(np.shape(x))
    maxi=np.empty([0])
    inds=np.empty([0])
    
    #Find the local maxima and the corresponding indices
    for ii in range(lx):
        #Pick candidate for maximum
        cand=x[ii]
        #Check that the candidate is above the threshold
        if cand<=lim:
            continue
        #Check whether maximum is a maximum cnt steps onward
        for jj in range(1,cnt+1):
            #Stop iteration if end of array is reached and
            #add the last found maximum, even if it's not a maximum
            #cnt steps onward, and is not eqautl to the previous maximum
            if jj+ii==lx:
                if cand>maxi[-1]:
                    maxi=np.concatenate([maxi,cand*np.ones([1])])
                    inds=np.concatenate([inds,ii*np.ones([1])])
                return maxi,inds
            #If a larger value is found
            #move the iteration to the new maximum value
            if x[ii+jj]>cand:
                ii=ii+jj
                break
            #If the same value is encountered again, discard it
            if len(maxi)>0:
                if cand==maxi[-1]:
                    continue
            #Save the maximum value and its index
            #and move the iteration to the next unchecked value
            if jj==cnt:    
                maxi=np.concatenate([maxi,cand*np.ones([1])])
                inds=np.concatenate([inds,ii*np.ones([1])])
                ii=ii+jj
    #Add the starting point of the series
    #maxi=np.concatenate([x[ii]*np.ones([1]),maxi])
    #inds=np.concatenate([ii*np.zeros([1]),inds])
    
    #Calculate and show autocorrelation of the declustered data points
    # if len(maxi)==0:
    #     autoc='nan'
    #     print("Length of maxis is 0.")
    #     print("Autocorrelation cannot be calculated")
    #     #return
    # else:
    #     autoc=sm.tsa.acf(maxi,nlags=cnt)
    # print(f"For declustering limit of {cnt} h the autocorrelation \
    #        with lag 1 is equal to {autoc}")
    
    return maxi,inds

#Calculate the standard deviation of given block size to
#define the erros for the BB method
def calc_errs(flux,d,ns,move_sd=None):
    #Convert d from days to hours
    d=int(d*24)
    
    errs=np.empty([0])
    lf=np.amax(np.shape(flux))
    num=int(np.floor(lf/d))
    
    
    if move_sd is None:
        for ii in range(num):
            subs=flux[ii*d:(ii+1)*d]
            var=np.var(subs)
            sd=var**0.5
            errs=np.concatenate((errs,ns*np.ones([d])*sd))
        
        #Add the final bit
        subs=flux[num*d:]
        var=np.var(subs)
        sd=var**0.5
        errs=np.concatenate((errs,ns*np.ones([len(subs)])*sd)) 
    else:
    #Calculate eror as a "moving s.d.", using the current point and the
    #preceding points, except in the beginning, where the s.d. is set constant
    #for the first d/2 data points If d is even, the last d/2 and the next d/2-1
    #values are included inthe calculation
        #subs=flux[lf-d:]
        d_flr=int(np.floor(d/2))
        subs=flux[:d_flr]
        var=np.var(subs)
        sd=var**0.5
        errs=np.concatenate((errs,ns*np.ones([len(subs)])*sd)) 
        # for ii in range(lf-d):
        #     subs=flux[ii:ii+d]
        #     var=np.var(subs)
        #     sd=var**0.5
        #     errs=np.concatenate((errs,np.ones([1])*sd))
        
        #Do errors for the data points from [:d/2] onwards, taking
        #the nearest d values around the data point
        for ii in range(lf-d):
            subs=flux[ii:ii+d]
            var=np.var(subs)
            sd=var**0.5
            errs=np.concatenate((errs,ns*np.ones([1])*sd))
        #Do the erros for the final d/2 data points
        subs=flux[lf-d_flr:]
        var=np.var(subs)
        sd=var**0.5
        errs=np.concatenate((errs,ns*np.ones([len(subs)])*sd)) 
        
        # except:
        #     subs=flux[lf-(d_ceil-1):]
        #     var=np.var(subs)
        #     sd=var**0.5
        #     errs=np.concatenate((errs,np.ones([len(subs)])*sd)) 
    return errs

#Calculate the iqr of a given block size to
#define the erros for the BB method
def calc_errs_iqr(flux,d,iq,move_sd=None):
    #Convert d from days to hours
    d=int(d*24)
    
    errs=np.empty([0])
    lf=np.amax(np.shape(flux))
    num=int(np.floor(lf/d))
    
    if move_sd is None:
        for ii in range(num):
            subs=flux[ii*d:(ii+1)*d]
            err=np.quantile(subs,iq)
            errs=np.concatenate((errs,err*np.ones([d])))
        
        #Add the final bit
        subs=flux[num*d:]
        err=np.quantile(subs,iq)
        errs=np.concatenate((errs,err*np.ones([len(subs)])))
    else:
    #Calculate error as a "moving iqr", using the current point and the
    #preceding points, except in the beginning, where the s.d. is set constant
    #for the first d/2 data points if d is even, the last d/2 and the next d/2-1
    #values are included inthe calculation
        #subs=flux[lf-d:]
        d_flr=int(np.floor(d/2))
        subs=flux[:d_flr]
        err=np.quantile(subs,iq)
        errs=np.concatenate((errs,err*np.ones([len(subs)])))
        # for ii in range(lf-d):
        #     subs=flux[ii:ii+d]
        #     var=np.var(subs)
        #     sd=var**0.5
        #     errs=np.concatenate((errs,np.ones([1])*sd))
        
        #Do errors for the data points from [:d/2] onwards, taking
        #the nearest d values around the data point
        for ii in range(lf-d):
            subs=flux[ii:ii+d]
            err=np.quantile(subs,iq)
            errs=np.concatenate((errs,err*np.ones([1])))
        #Do the erros for the final d/2 data points
        subs=flux[lf-d_flr:]
        err=np.quantile(subs,iq)
        errs=np.concatenate((errs,err*np.ones([len(subs)])))
        
    return errs
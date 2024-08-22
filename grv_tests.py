#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 13:50:23 2021

@author: mikkosavola
"""
#Script for testing the functions in file fat_tail_stats.py
#using GRVs
import aux_mi as am
import main_mi as mm
import numpy as np
import scipy.stats as stats
import os
import numpy.random as nr
import fat_tail_stats as fts

#Sample size is N
N=10**4
noise=10**(-4)
rans=fts.add_noise(nr.normal(size=N),noise)
#Set k Pickand's and Hill plots
k=400
#Set the maximum u for m2s plots
m2s_n_x=100
#Define a new standard deviation but keep the distribution centralised
sigma=10
#Move lower boundary to zero.
rans=rans*sigma
mini=abs(min(rans))
rans=rans+mini

m2s_p=[1,2,3,4]
ttl='GRV_'+'N='+str(N)+'_'+'sigma='+str(sigma)

hdr='Sample values from N(0,1) for n='+str(N)
fts.do_plots(rans,k,ttl,hdr,m2s_p,m2s_n_x=m2s_n_x)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 14:38:31 2021

@author: mikkosavola
"""
#Script for testing the functions in file fat_tail_stats.py
#using a power law distribution
import aux_mi as am
import main_mi as mm
import numpy as np
import scipy.stats as stats
import os
import numpy.random as nr
import fat_tail_stats as fts
import matplotlib.pyplot as plt

#Shape and scale parameters of the power distributions
a1=1
a2=2
a3=3
m1=1
m2=1
m3=1
a4=0.5
m4=1
a5=0
m5=1
a6=-0.5
m6=1
#Support of the sample
#s0=0
#s1=10**4
#Sample size is N
N=10**4
#NUmber of runs for testing the average xi and its confidenve interval
n=10**3
conf=0.95

#Create samples from uniform distributions
#See https://mathworld.wolfram.com/RandomNumber.html
# unis1=nr.uniform(0,1,N)
# unis2=nr.uniform(0,1,N)
#unis3=nr.uniform(0,1,N)
# min1=np.amin(rans1)
# min2=np.amin(rans2)
# min3=np.amin(rans3)
# max1=np.amax(rans1)
# max2=np.amax(rans2)
# max3=np.amax(rans3)
#Define constants to amke the samplew distributions integrable to 1
# C1=(a1+1)/(s1**(a1+1)-s0**(a1+1))
# C2=(a2+1)/(s1**(a2+1)-s0**(a2+1))
# C3=(a3+1)/(s1**(a3+1)-s0**(a3+1))

#Create samples from classical Pareto distributions
#by adding one to the Lomax dist. and multiplying by m
# rng=nr.default_rng()
# rans1=(rng.pareto(a1,N)+1)*m1
# rans2=(rng.pareto(a2,N)+1)*m2
# rans3=(rng.pareto(a3,N)+1)*m3

#Use scipys generator for GPD, instead of the above classical Pareto
GPD=lambda a,scale,size : stats.genpareto.rvs(a,scale=scale,size=size)
rans1=GPD(a1,scale=m1,size=N)
rans2=GPD(a2,scale=m2,size=N)
rans3=GPD(a3,scale=m3,size=N)
#rans1=np.reshape(rans1,[N,1])
#rans2=np.reshape(rans2,[N,1])
#rans3=np.reshape(rans3,[N,1])

#Set k for Pickand's and Hill plots
k=400
if 4*k>N:
    k=np.floor(N/4)-1
    k=int(k)
    
#Set the maximum u for m2s plots
m2s_n_x=100

m2s_p=[1,2,3,4]
ttl1='GPD with shape a='+str(a1)+' and scale m='+str(m1)
ttl2='GPD with shape a='+str(a2)+' and scale m='+str(m2)
ttl3='GPD with shape a='+str(a3)+' and scale m='+str(m3)
ttl4='GPD with shape a='+str(a4)+' and scale m='+str(m4)
ttl5='GPD with shape a='+str(a5)+' and scale m='+str(m5)
ttl6='GPD with shape a='+str(a6)+' and scale m='+str(m6)

hdr1='Sample values from GPD(a='+str(a1)+', m='+str(m1)+')'
hdr2='Sample values from GPD(a='+str(a2)+', m='+str(m2)+')'
hdr3='Sample values from GPD(a='+str(a3)+', m='+str(m3)+')'
fts.do_plots(rans1,k,ttl1,hdr1,m2s_p,m2s_n_x=m2s_n_x)
fts.do_plots(rans2,k,ttl2,hdr2,m2s_p,m2s_n_x=m2s_n_x)
fts.do_plots(rans3,k,ttl3,hdr3,m2s_p,m2s_n_x=m2s_n_x)
# fts.pick_plot_N(GPD,a1,m1,N,k,ttl1,n,conf)
# fts.hill_plot_N(GPD,a1,m1,N,k,ttl1,n,conf)
# fts.pick_plot_N(GPD,a2,m2,N,k,ttl2,n,conf)
# fts.hill_plot_N(GPD,a2,m2,N,k,ttl2,n,conf)
# fts.pick_plot_N(GPD,a3,m3,N,k,ttl3,n,conf)
# fts.hill_plot_N(GPD,a3,m3,N,k,ttl3,n,conf)
# fts.pick_plot_N(GPD,a4,m4,N,k,ttl4,n,conf)
# fts.hill_plot_N(GPD,a4,m4,N,k,ttl4,n,conf)
# fts.pick_plot_N(GPD,a5,m5,N,k,ttl5,n,conf)
# fts.hill_plot_N(GPD,a5,m5,N,k,ttl5,n,conf)
# fts.pick_plot_N(GPD,a6,m6,N,k,ttl6,n,conf)
# fts.hill_plot_N(GPD,a6,m6,N,k,ttl6,n,conf)

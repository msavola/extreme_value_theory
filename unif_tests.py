#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 14:30:49 2021

@author: mikkosavola
"""
#Script for testing the functions in file fat_tail_stats.py
#using Uniform distributions
import aux_mi as am
import main_mi as mm
import numpy as np
import scipy.stats as stats
import os
import numpy.random as nr
import fat_tail_stats as fts
import matplotlib.pyplot as plt

#Limits of the uniform distributions
a1=0
b1=1
a2=1
b2=100
a3=14
b3=27
#Sample size is N
N=10**4
#Create samples from the uniform distribution
rans1=nr.uniform(a1,b1,N)
rans2=nr.uniform(a2,b2,N)
rans3=nr.uniform(a3,b3,N)
#Set k Pickand's and Hill plots
k=400
if 4*k>N:
    k=np.floor(N/4)-1
    k=int(k)
    
#Set the maximum u for m2s plots
m2s_n_x=100


m2s_p=[1,2,3,4]
ttl1='Unif('+str(a1)+','+str(b1)+')_'+'N='+str(N)+')'
ttl2='Unif('+str(a2)+','+str(b2)+')_'+'N='+str(N)+')'
ttl3='Unif('+str(a3)+','+str(b3)+')_'+'N='+str(N)+')'

hdr1='Sample values from Uniform('+str(a1)+','+str(b1)+') for n='+str(N)
hdr2='Sample values from Uniform('+str(a2)+','+str(b2)+') for n='+str(N)
hdr3='Sample values from Uniform('+str(a3)+','+str(b3)+') for n='+str(N)
fts.do_plots(rans1,k,ttl1,hdr1,m2s_p,m2s_n_x=m2s_n_x)
fts.do_plots(rans2,k,ttl2,hdr2,m2s_p,m2s_n_x=m2s_n_x)
fts.do_plots(rans3,k,ttl3,hdr3,m2s_p,m2s_n_x=m2s_n_x)



#Set integer values for u
us1=np.linspace(a1,m2s_n_x,m2s_n_x-a1+1)
us2=np.linspace(a2,m2s_n_x,m2s_n_x-a2+1)
us3=np.linspace(a3,m2s_n_x,m2s_n_x-a3+1)
#Evaluate the mean excess function analytically for different
#uniform distributions
e1=0.5*(b1-us1)
e2=0.5*(b2-us2)
e3=0.5*(b3-us3)
#Plot analytical me plots
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.plot(us1,e1,label=ttl1)
ax.plot(us2,e2,label=ttl2)
ax.plot(us3,e3,label=ttl3)
ax.set_xlabel('threshold')
ax.set_ylabel('mean excess')
ax.set_title('Analytical mef of uniform distributions')
fig.legend()
fig.show()
fig.savefig('me_unif_analytical_test.pdf')
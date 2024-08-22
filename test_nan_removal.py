#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 15:41:59 2021

@author: mikkosavola
"""

import numpy as np
import fat_tail_stats as fts

test1=np.array([float('nan'),1,2,3,4])
test2=np.array([1991,1991,float('nan'),1993,1994])
big=9999
new1,new2=fts.clean_data(test1,test2,big)


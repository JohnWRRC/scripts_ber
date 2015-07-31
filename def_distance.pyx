# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 10:46:39 2015

@author: John
"""

import numpy as np

def distance(xy_ind_a, xy_ind_b):
    
    cdef float a_x, a_y, b_x ,b_y,dist
    
    a_x=xy_ind_a[0]
    a_y=xy_ind_a[1]
    b_x=xy_ind_b[0]
    b_y=xy_ind_b[1]

    dist = np.sqrt((a_x-b_x)*(a_x-b_x) + (a_y-b_y)*(a_y-b_y))
    return dist
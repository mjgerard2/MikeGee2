#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:53:47 2019

@author: michael
"""

import numpy as np

def line_val(p1, p2, x):
    x1, y1 = p1
    x2, y2 = p2
    
    m = (y2 - y1) / (x2 - x1)
    return m*x - m*x1 + y1
    
p1 = [.4, .09276]
p2 = [.5, .08445]
x = .4118
print(line_val(p1, p2, x))
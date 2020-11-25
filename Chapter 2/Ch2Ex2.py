# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 19:18:17 2020

@author: LabUser
"""



import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
import MeriodionalPlaneOptics as meri

r1 = 10
r2 = -10
t1 = 2
n1 = 1
n2 = 1.5
lp = -1200
alpha1 = 0
c1 = r1-lp
x = np.arange(0.00001,4.1,0.1)

lpdash = np.zeros(x.shape)
ALSA = np.zeros(x.shape)
TSA = np.zeros(x.shape)
for ii, entry in enumerate(x):
    T1 = meri.rayDisplacement(alpha1, entry, c1, r1)

    K1 = meri.skewPower(alpha1, entry, c1, r1, n1, n2)

    T1P = meri.TranslationMatrix(T1, n1)

    R1 = meri.RefractionMatrix(K1)
    
    result1 = R1.dot(T1P.dot(np.array([[0],[entry]])))
    
    x1 = result1[1,0]
    
    sin_alpha2 = result1[0,0]/n2
    alpha2 = math.asin(sin_alpha2)
    
    
    #repeat everything with shift in origin
    c2 = meri.shiftOrigin(T1, alpha1, c1, r1, r2, t1)
    T2 = meri.rayDisplacement(alpha2,x1,c2,r2)
    K2 = meri.skewPower(alpha2, x1, c2, r2, n2, n1)
    T21 = meri.TranslationMatrix(T2, n2)
    R2 = meri.RefractionMatrix(K2)
    result2 = R2.dot(T21.dot(np.array([[n2*math.sin(alpha2)],[x1]])))
    sinalpha2dash = result2[0,0]/n1
    alpha2dash = math.asin(sinalpha2dash)
    x2 = result2[1,0]
    lpdash[ii] = x2/(math.tan(-alpha2dash))-(-r2-math.sqrt(r2**2-x2**2))
    ALSA[ii] = lpdash[ii]-lpdash[0]
    TSA[ii] = ALSA[ii]*math.tan(-alpha2dash)
print(lpdash)
print(ALSA)
print(TSA)
print(x)
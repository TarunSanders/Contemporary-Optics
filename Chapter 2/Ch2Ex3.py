# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:03:16 2020

@author: LabUser
"""


import numpy as np
import MeriodionalPlaneOptics as meri
import math
import matplotlib.pyplot as plt

r1 = np.array([-10, 1000000,20,10,6.67,5,3.33])
r2 = np.array([-3.33,-5,-6.67,-10,-20,1000000,10])
S = np.array([-2,-1,-0.5,0,0.5,1,2])
ALSA = np.zeros(S.shape)
t=2
n = 1.5
lp=-1200
x = np.array([0.00001,1])
lpdash = np.zeros(x.shape)

j = 0
for radius1,radius2 in zip(r1,r2):
    i = 0
    
    for height in x:
        #ray tracing to surface 1
        c1 = radius1 - lp
        T = meri.rayDisplacement(0,height,c1,radius1)
        K = meri.skewPower(0,height,c1,radius1,1,n)
        T1P = meri.TranslationMatrix(T,1)
        R1 = meri.RefractionMatrix(K)
        S1P = R1.dot(T1P)
        result1 = S1P.dot(np.array([[1*math.sin(0)],[height]]))
        x1dash = result1[1,0] 
        alpha1 = math.asin(result1[0,0]/n)
        
        
        #ray tracing to surface 2
        c2 = meri.shiftOrigin(T,0,c1,radius1,radius2,t)
        T = meri.rayDisplacement(alpha1,x1dash,c2,radius2)
        K = meri.skewPower(alpha1,x1dash,c2,radius2,n,1)
        T21= meri.TranslationMatrix(T,n)
        R2 = meri.RefractionMatrix(K)
        S21 = R2.dot(T21)
        
        result2 = S21.dot(np.array([[n*math.sin(alpha1)],[x1dash]]))
        x2 = result2[1,0]
        alpha2 = math.asin(result2[0,0])
        lpdash[i] = x2 / math.tan(-alpha2) - (abs(radius2)-math.sqrt(radius2**2-x2**2))
        i += 1
    ALSA[j] = lpdash[1]-lpdash[0]
    j += 1

plt.plot(S,ALSA,'o',color='black')
plt.xlabel("Shape factor, S",fontsize=20)
plt.ylabel("Longitudinal aberration, ALSA",fontsize=20)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.xlim(-2,2)
plt.ylim(0,-1.2)
plt.show()

r = np.arange(0,1.05,0.001)
ALSAnorm = np.array([(entry**2)*(entry**2-1) for entry in r])
plt.plot(ALSAnorm,r)

plt.ylim(0,1.1)
plt.xlim(-0.3,0.2)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.xlabel("Normalised ALSA, ALSA'",fontsize=20)
plt.ylabel("Normalised object height, r",fontsize=20)
plt.show()
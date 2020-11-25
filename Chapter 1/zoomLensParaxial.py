# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 17:24:04 2020

@author: LabUser
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

def systemMatrix(t,x):
    S21 = np.array([[1,-1/10],[0,1]])
    T32 = np.array([[1,0],[t,1]])
    S43 = np.array([[1,-1/20],[0,1]])
    T54 =np.array([[1,0],[60-t+x,1]])
    S65 = np.array([[1,-1/10],[0,1]])
    S61 = S65.dot(T54.dot(S43.dot((T32.dot(S21)))))
    return S61

def objectToImageDistance(S61,x):
    Limage =(-100*S61[1,1]-S61[1,0])/(100*S61[0,1]+S61[0,0])
    D = 160 + x + Limage              
    return D

def magnification(S61):
    Beta = 1/(S61[0,0]+100*S61[0,1])
    return Beta

def focalLength(S61):
    fdash = -1/S61[0,1]
    return fdash

t = np.linspace(0,60,100)
sysMatrices = [systemMatrix(entry,0) for entry in t]

Beta = np.array([magnification(entry) for entry in sysMatrices])

D =  np.array([objectToImageDistance(entry,0) for entry in sysMatrices])

fdash = np.array([focalLength(entry) for entry in sysMatrices])



plt.plot(t,D)
plt.xlim(0,60)
plt.ylim(150,180)
plt.xlabel("Separation of lens 1 and lens 2, t",fontsize=20)
plt.ylabel("Object-to-image distance, D",fontsize=20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.show()

plt.plot(t,Beta)
plt.xlim(0,60)
plt.ylim(0,0.06)
plt.xlabel("Separation of lens 1 and lens 2, t",fontsize=20)
plt.ylabel("Magnification",fontsize=20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.show()

plt.plot(t,fdash)
plt.show()

#movement of lens 3 to correct for defocusing
S61 = systemMatrix(30,0)
Dtrue = objectToImageDistance(S61,0)
t =  np.linspace(30,60,100)
error = 0.1
spacing = np.zeros(len(t))
for i, entry in enumerate(t):
    x = 0
    S61 = systemMatrix(entry,x)
    D = objectToImageDistance(S61,x)
    while abs(Dtrue-D) > 0.1:
        x+=0.001
        S61 = systemMatrix(entry,x)
        D = objectToImageDistance(S61,x)
    spacing[i] = 60 + x
    

plt.plot(t,spacing) 
plt.xlim(30,60)
plt.ylim(60,63)
plt.xlabel("Separation of lens 1 and lens 2, t",fontsize=15)
plt.ylabel("Separation between lens 1 and lens 3",fontsize=15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.text(40,61,"D = 170.3 (+- 0.1)", fontsize = 15)
plt.show()
    

    


# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 15:23:14 2020

@author: LabUser
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

def rayDisplacement(alpha,x,c, r):
    if r > 0:
        T = c*math.cos(alpha)-x*math.sin(alpha)-math.sqrt(r**2-(x*math.cos(alpha)+c*math.sin(alpha))**2)
    elif r < 0:
        T = c*math.cos(alpha)-x*math.sin(alpha)+math.sqrt(r**2-(x*math.cos(alpha)+c*math.sin(alpha))**2)
    else:
        print("Invalid curvature radius of surface")
    return T

def skewPower(alpha,x,c,r,n1,n2):
    sin_Theta = (x*math.cos(alpha)+c*math.sin(alpha))/r
    
    Theta = math.asin(sin_Theta)
    K = (math.sqrt(n2**2-(n1**2)*(sin_Theta**2))-n1*math.cos(Theta))/r
    return K

def TranslationMatrix(T,n):
     T21 = np.array([[1,0],[T/n, 1]])
     return T21
 
def RefractionMatrix(K):
    R = np.array([[1,-K],[0, 1]])
    return R
    
def shiftOrigin(Told,alphaOld,cold,r1,r2,t1):
    cnew = -Told*math.cos(alphaOld)+cold-r1+t1+r2
    return cnew

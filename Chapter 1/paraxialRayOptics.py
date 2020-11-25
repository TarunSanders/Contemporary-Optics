# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 17:53:42 2020

@author: LabUser
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

#left is object side and right is image side in these definitions
def Refraction(nleft,nright,r):
    k = (nright-nleft)/r
    R = np.array([[1,-k],[0, 1]])
    return R

def Translation(t,n):
    T = np.array([[1,0],[t/n, 1]])
    return T

def SystemMatrix(Rright,T21,Rleft):
    S21 = Rright.dot(T21.dot(Rleft))
    return S21

def cardinalPlanes(S,nleft,nright):
  # Gaussian constants from system matrix
    a = -S[0,1]
    b = S[0,0]
    c = S[1,1]
    d = -S[1,0]
    # r refers to image side, H for unit planes, l for distance from vertices (lens edge), F for focal plane, f for focal length (unit to focal plane distance)
    lHr = nright*(c-1)/a
    lH = nleft*(1-b)/a
    lF = -nleft*b/a
    lFr = nright*c/a
    f = lF-lH
    fr = lFr-lHr
    return lH,lHr,lF,lFr, f, fr

#this function returns image position lp' (distance from last lens) and objectTo image distance and magnification
def image(objectPosition,S,nleft,nright):
    c = S[1,1]
    b = S[0,0]
    a = -S[0,1]
    M = 1/(b+a*objectPosition/nleft)
    imagePosition = nright*(c-M)/a
    objectToimageDistance = imagePosition - objectPosition
    return imagePosition, objectToimageDistance, M
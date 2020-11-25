# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:21:55 2020

@author: LabUser
"""

import numpy as np
import SkewRayOptics as skew
import math
import matplotlib.pyplot as plt


#lens and medium parameters
r1 = 10
r2 = -10
n1 = 1
n1dash = 1.5
t = 2

#definition of rays
cx1top = 0.01
cy1top = 0
cz1top = math.sqrt(1-cx1top**2-cy1top**2)

x1top = 1
y1top = 0
z1top = r1 - math.sqrt(r1**2 - x1top**2)

cx1bot = cx1top
cy1bot = 0
cz1bot = math.sqrt(1-cx1top**2-cy1top**2)

x1bot = -x1top
y1bot = -y1top
z1bot = r1 - math.sqrt(r1**2 - x1bot**2)




def rayLens(x1,y1,z1,cx1,cy1,cz1,n1,n1dash,t,r1,r2):
    lpdash = z1
    cosTheta1 = skew.cosTheta(x1, y1, z1, lpdash, r1, cx1, cy1, cz1)
    K1 = skew.skewPower(n1, n1dash, r1, cosTheta1)
    R1 = skew.RefractionMatrix(K1)
    result1 = R1.dot(np.array([[n1*cx1,n1*cy1,n1*cz1],[x1,y1,z1]]))
    x1dash,y1dash,z1dash = result1[1]
    cx1dash = result1[0,0]/n1dash
    cy1dash = result1[0,1]/n1dash
    cz1dash = (result1[0,2]+K1*r1)/n1dash

    lpdash = z1dash - t #left of second surface
    #translation
    cosTheta2 = skew.cosTheta(x1dash, y1dash, z1dash, lpdash, r2, cx1dash, cy1dash, cz1dash)
    T2 = skew.rayDisplacement(x1dash, y1dash, z1dash, lpdash, r2, cx1dash, cy1dash, cz1dash, cosTheta2)
    K2 = skew.skewPower(n1dash, n1, r2, cosTheta2)
    T21 = skew.TranslationMatrix(n1dash, T2)
    trans = T21.dot(np.array([[n1dash*cx1dash,n1dash*cy1dash,n1dash*cz1dash],[x1dash,y1dash,z1dash]]))
    x2, y2, z2 = trans[1]
    z2 = z2 - t
    R2 = skew.RefractionMatrix(K2)
    result2 = R2.dot(np.array([[n1dash*cx1dash,n1dash*cy1dash,n1dash*cz1dash],[x2,y2,z2]]))
    x2dash,y2dash,z2dash = result2[1]
    cx2dash = result2[0,0] / n1
    cy2dash = result2[0,1] / n1
    cz2dash = (result2[0,2]+K2*r2) / n1
    return x2dash,y2dash,z2dash, cx2dash, cy2dash,cz2dash

#basically solving two simultaneous linear equations see Ch. 3 doc 
def findIntersection(x2dashtop,x2dashbot,z2dashtop,z2dashbot,cx2dashtop,cx2dashbot,cz2dashtop,cz2dashbot):
    mtop = cx2dashtop / cz2dashtop #gradient of top ray
    mbot = cx2dashbot / cz2dashbot #gradient of bottom ray
    ctop = x2dashtop - mtop*z2dashtop # interseept with plane with V2
    cbot = x2dashbot - mbot*z2dashbot
    FP = (cbot - ctop)/(mtop-mbot)
    return FP

def rayFPplane(FP,x2dash,y2dash,z2dash, cx2dash, cy2dash,cz2dash,n1):
    lpdash = FP - z2dash
    T = lpdash / cz2dash
    TFP2 = skew.TranslationMatrix(n1, T)
    trans = TFP2.dot(np.array([[n1*cx2dash,n1*cy2dash,n1*cz2dash],[x2dash,y2dash,z2dash]]))
    x3, y3, z3 = trans[1]
    z3 = z3 - FP
    return x3,y3,z3

#emerging top ray parameters
x2dashtop,y2dashtop,z2dashtop, cx2dashtop, cy2dashtop,cz2dashtop = rayLens(x1top,y1top,z1top,cx1top,cy1top,cz1top,n1,n1dash,t,r1,r2)
#emerging bottom ray parameters
x2dashbot,y2dashbot,z2dashbot, cx2dashbot, cy2dashbot,cz2dashbot = rayLens(x1bot,y1bot,z1bot,cx1bot,cy1bot,cz1bot,n1,n1dash,t,r1,r2)

FP = findIntersection(x2dashtop, x2dashbot, z2dashtop, z2dashbot, cx2dashtop, cx2dashbot, cz2dashtop, cz2dashbot)

xFP,yFP,zFP = rayFPplane(FP, x2dashtop, y2dashtop, z2dashtop, cx2dashtop, cy2dashtop, cz2dashtop, n1)
print(FP)
print(np.array([xFP,yFP,zFP]))




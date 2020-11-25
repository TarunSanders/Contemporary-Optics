# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 09:52:47 2020

@author: LabUser
"""

import numpy as np
import SkewRayOptics as skew
import math
import matplotlib.pyplot as plt

r1 = 50
r2 = -50
x = 10
y = 1
z = 2

cx = -0.1
cy = 0.1
cz = math.sqrt(1-cx**2-cy**2)

lpdash = -8

n1 = 1
n1Dash = 1.5

t1 = 15 #lens thickness

FP = 65.517241 # paraxial focus to the right of lens #use prev meriodonal plane or paraxial approx code to determine
r3 = math.inf 
# first surface translation and refraction
cosTheta1 = skew.cosTheta(x, y, z, lpdash, r1, cx, cy, cz)
T1 = skew.rayDisplacement(x, y, z, lpdash, r1, cx, cy, cz, cosTheta1)
T1P = skew.TranslationMatrix(n1, T1)
K1 = skew.skewPower(n1, n1Dash, r1, cosTheta1)
R1 = skew.RefractionMatrix(K1)
    #super matrix multiplication for x,y, and z
Translation1 = T1P.dot(np.array([[n1*cx, n1*cy, n1*cz],[x, y, z]]))
Translation1[1,2] = Translation1[1,2] - (z - lpdash)
result1 = R1.dot(Translation1)

cx1Dash = result1[0,0] / n1Dash
cy1Dash = result1[0,1] / n1Dash
cz1Dash = (result1[0,2] + K1*r1)/ n1Dash

x1Dash = result1[1,0] 
y1Dash = result1[1,1] 
z1Dash = result1[1,2] 

# second surface translation and refraction
cosTheta2 = skew.cosTheta(x1Dash, y1Dash, z1Dash, -(t1 - z1Dash), r2, cx1Dash, cy1Dash, cz1Dash)
T2 = skew.rayDisplacement(x1Dash, y1Dash, z1Dash, -(t1 - z1Dash), r2, cx1Dash, cy1Dash, cz1Dash, cosTheta2)
T21 = skew.TranslationMatrix(n1Dash, T2)
K2 = skew.skewPower(n1Dash, n1, r2, cosTheta2)
R2 = skew.RefractionMatrix(K2)

Translation2 = T21.dot(np.array([[n1Dash*cx1Dash, n1Dash*cy1Dash, n1Dash*cz1Dash],[x1Dash, y1Dash, z1Dash]]))
Translation2[1,2] = Translation2[1,2] - t1
result2 = R2.dot(Translation2)

cx2Dash = result2[0,0] / n1
cy2Dash = result2[0,1] / n1
cz2Dash = (result2[0,2] + K2*r2)/ n1

x2Dash = result2[1,0] 
y2Dash = result2[1,1] 
z2Dash = result2[1,2] 

#position of ray at paraxial focal plane
T3 = (FP-z2Dash)/cz2Dash
T32 = skew.TranslationMatrix(n1, T3)
result3 = T32.dot(np.array([[n1*cx2Dash, n1*cy2Dash, n1*cz2Dash],[x2Dash, y2Dash, z2Dash]]))
x3 = result3[1,0]
y3Dash = result3[1,1] 
z3Dash = result3[1,2] - FP
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 15:42:23 2020

@author: LabUser
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
import MeriodionalPlaneOptics as meri

r1 = 50
r2 = -50
t1 = 15
n2 = 1.5
n1 = 1
lp = -200
sin_alpha = 0.1
x = 0
c1 = r1-lp

alpha = math.asin(sin_alpha)
print(c1)
#first ray displacement from object to surface
T1 = meri.rayDisplacement(alpha, x, c1, r1)
#first skew power at surface
K1 = meri.skewPower(alpha, x, c1, r1, n1, n2)

print("T and K are:")
print(T1)
print(K1)

#First translation matrix
T1P = meri.TranslationMatrix(T1, n1)
#First refraction matrix
R1 = meri.RefractionMatrix(K1)

result1 = R1.dot(T1P.dot(np.array([[n1*sin_alpha],[x]])))


x1 = result1[1,0]
sin_alpha2 = result1[0,0]/n2
alpha2 = math.asin(sin_alpha2)

print(x1)
print(n2*math.sin(alpha2))

#repeat everything with shift in origin
c2 = meri.shiftOrigin(T1, alpha, c1, r1, r2, t1)
print("c2 is")
print(c2)
T2 = meri.rayDisplacement(alpha2,x1,c2,r2)

K2 = meri.skewPower(alpha2, x1, c2, r2, n2, n1)
print("K")
print(K2)
T21 = meri.TranslationMatrix(T2, n2)
R2 = meri.RefractionMatrix(K2)
result2 = R2.dot(T21.dot(np.array([[n2*math.sin(alpha2)],[x1]])))

#x2 = x1+T2*math.sin(alpha2)
x2 = result2[1,0]
sin_alpha2dash = result2[0,0]/n1
#sin_alpha2dash = (n2*math.sin(alpha2)-K2*x2)/n1
alpha2dash = math.asin(sin_alpha2dash)
print(x2)
print(alpha2dash)

#finding location wher ray crosses axis
lpdash = x2/(math.tan(-alpha2dash))-(-r2-math.sqrt(r2**2-x2**2))
print(lpdash)
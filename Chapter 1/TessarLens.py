# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 14:22:27 2020

@author: LabUser
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
import paraxialRayOptics as para

#computing first system matrix (aysymmetric double convex lens)
R1 = para.Refraction(1, 1.6116, 1.628)
# print(R1)
T21 = para.Translation(0.357, 1.6116)
R2 = para.Refraction(1.6116,1,-27.57)
#print(R2)
S21 = para.SystemMatrix(R2,T21,R1) 
#print(S21)

#Translation in air 1
T32 = para.Translation(0.189, 1)

#computing element system matrix (asymetric double concave lens)
R3 = para.Refraction(1, 1.6053, -3.457)
T43 = para.Translation(0.081,1.6053)
R4 = para.Refraction(1.6053,1,1.582)
S43 = para.SystemMatrix(R4, T43, R3)

#Translation in air 2
T54 = para.Translation(0.325,1)

#system matrix for planoconcave part of doublet
R5 = para.Refraction(1,1.5123,float('inf'))
T65 = para.Translation(0.217, 1.5123)
R6 = para.Refraction(1.5123, 1.6116, 1.920)
S65 = para.SystemMatrix(R6, T65, R5)

#system matrix for asymmettric convex part of doublet
T76 = para.Translation(0.396, 1.6116)
R7 = para.Refraction(1.6116,1,-2.4)
S76 = R7.dot(T76)

#print(S76)
#entire system matrix
S71 = S76.dot(S65.dot(T54.dot(S43.dot(T32.dot(S21)))))
lH,lHr,lF,lFr, f, fr = para.cardinalPlanes(S71, nleft=1, nright=1) #air on both sides of system
print("lH,lHr,lF,lFr, f, fr are in order: ")
print(lH)
print(lHr)
print(lF)
print(lFr)
print(f)
print(fr)

lPr, D, M = para.image(-100,S71,1,1)

print("image position to right of last lens, object-to-image distance and magnification are:")
print(lPr)
print(D)
print(M)
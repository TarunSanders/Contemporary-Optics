# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:28:19 2020

@author: LabUser
"""



import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
import paraxialRayOptics as para

R1 = para.Refraction(1,1.5,10)
T21 = para.Translation(2,1.5)
R2 = para.Refraction(1.5,1,-10)

S21 = para.SystemMatrix(R2, T21, R1)
lH,lHr,lF,lFr, f, fr = para.cardinalPlanes(S21, 1,1)

print("Paraxial focal plane location is:")
print(lFr)
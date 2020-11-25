# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:06:52 2020

@author: LabUser
"""


import numpy as np
import SkewRayOptics as skew
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator)

#lens and medium parameters
r1 = 10
r2 = -10
n1 = 1
n1dash = 1.5
t = 2


#finds position relative to V2 and directional cosines of emerging ray from lens
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

x1 = np.arange(0.1,1.1,0.1)
y1 = 0
z1 = r1-np.sqrt(r1**2-x1**2)
cx1 = np.arange(0.001,0.011,0.001)
cy1 = 0
cz1 = np.sqrt(1-cx1**2)
FP = np.zeros((x1.shape[0],cx1.shape[0]))

#finding matrix of focal planes relative to v2 for 100 combinations of angles and zone radius
for i in range(x1.shape[0]):
    for j in range(cx1.shape[0]):
        x2dashtop,y2dashtop,z2dashtop,cx2dashtop, cy2dashtop,cz2dashtop = rayLens(x1[i], y1, z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        x2dashbot,y2dashbot,z2dashbot,cx2dashbot, cy2dashbot,cz2dashbot = rayLens(-x1[i], y1, z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        FP[i,j] = findIntersection(x2dashtop, x2dashbot, z2dashtop, z2dashbot, cx2dashtop, cx2dashbot, cz2dashtop, cz2dashbot)
print(FP)

r = x1 #zone radius
#phi = math.pi/180*np.arange(0,181,15)
A = np.zeros((x1.shape[0],cx1.shape[0]))
Cs = np.zeros((x1.shape[0],cx1.shape[0]))
Ct = np.zeros((x1.shape[0],cx1.shape[0]))
Xmax = np.zeros((x1.shape[0],cx1.shape[0]))
Xmin = np.zeros((x1.shape[0],cx1.shape[0]))
XL = np.zeros(cx1.shape)
#X = np.zeros(phi.shape)
#Y = np.zeros(phi.shape)
S = np.zeros((x1.shape[0],cx1.shape[0]))

for i in range(r.shape[0]):
    for j in range(cx1.shape[0]):
        # for k in range(phi.shape[0]):
        #     x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(r[i]*math.sin(phi[k]), r[i]*math.cos(phi[k]), z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        #     x3,y3,z3 = rayFPplane(FP[i,j], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)
        #     X[k] = x3
        #max x occurs at x=zone radius and y=0
        x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(r[i], 0, z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        x3,y3,z3 = rayFPplane(FP[i,j], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)
        Xmax[i,j] = x3*1000000 #amplified value
        #min x occurs ay x = 0, and y = zone radius
        x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(0, r[i], z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        x3,y3,z3 = rayFPplane(FP[i,j], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)
        Xmin[i,j] = x3*1000000
        
        if i == 0:
            XL[j] = Xmax[i,j]


for j in range(cx1.shape[0]):
    for i in range(r.shape[0]):
        Cs[i,j] = (Xmax[i,j] - Xmin[i,j])/2
        Ct[i,j] = Xmin[i,j] - XL[j]
        A[i,j] = Ct[i,j]*Cs[i,j]
        S[i,j] = Cs[i,j]/(Ct[i,j])
        print([Xmax[i,j],Xmin[i,j],XL[j],Ct[i,j],Cs[i,j],A[i,j],S[i,j]])
        
        
r1 = 6.67
r2 = -20
n1 = 1
n1dash = 1.5
t = 2

x1 = np.arange(0.1,1.1,0.1)
y1 = 0
z1 = r1-np.sqrt(r1**2-x1**2)
for i in range(x1.shape[0]):
    for j in range(cx1.shape[0]):
        x2dashtop,y2dashtop,z2dashtop,cx2dashtop, cy2dashtop,cz2dashtop = rayLens(x1[i], y1, z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        x2dashbot,y2dashbot,z2dashbot,cx2dashbot, cy2dashbot,cz2dashbot = rayLens(-x1[i], y1, z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        FP[i,j] = findIntersection(x2dashtop, x2dashbot, z2dashtop, z2dashbot, cx2dashtop, cx2dashbot, cz2dashtop, cz2dashbot)
        
for i in range(r.shape[0]):
    for j in range(cx1.shape[0]):
        # for k in range(phi.shape[0]):
        #     x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(r[i]*math.sin(phi[k]), r[i]*math.cos(phi[k]), z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        #     x3,y3,z3 = rayFPplane(FP[i,j], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)
        #     X[k] = x3
        #max x occurs at x=zone radius and y=0
        x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(r[i], 0, z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        x3,y3,z3 = rayFPplane(FP[i,j], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)
        Xmax[i,j] = x3*1000000 #amplified value
        #min x occurs ay x = 0, and y = zone radius
        x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(0, r[i], z1[i], cx1[j], cy1, cz1[j], n1, n1dash, t, r1, r2)
        x3,y3,z3 = rayFPplane(FP[i,j], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)
        Xmin[i,j] = x3*1000000
        
        if i == 0:
            XL[j] = Xmax[i,j]

print('Second lens with minimal spherical aberration:')
for j in range(cx1.shape[0]):
    for i in range(r.shape[0]):
        Cs[i,j] = (Xmax[i,j] - Xmin[i,j])/2
        Ct[i,j] = Xmin[i,j] - XL[j]
        A[i,j] = Ct[i,j]*Cs[i,j]
        S[i,j] = Cs[i,j]/(Ct[i,j])
        print([Xmax[i,j],Xmin[i,j],XL[j],Ct[i,j],Cs[i,j],A[i,j], S[i,j]])

phi = math.pi/180*np.arange(0,181,15)
XPHI = np.zeros((phi.shape[0],r.shape[0]))
YPHI = np.zeros((phi.shape[0],r.shape[0]))
#Z = np.zeros((phi.shape[0],r.shape[0]))
cx1 = 0.001
cy1 = 0
cz1 = np.sqrt(1-cx1**2)

for i in range(r.shape[0]):
    for j in range(phi.shape[0]):
        x2dash,y2dash,z2dash,cx2dash, cy2dash,cz2dash = rayLens(r[i]*math.sin(phi[j]), r[i]*math.cos(phi[j]), z1[i], cx1, cy1, cz1, n1, n1dash, t, r1, r2)
        XPHI[j,i],YPHI[j,i],z3 = rayFPplane(FP[i,0], x2dash, y2dash, z2dash, cx2dash, cy2dash, cz2dash, n1)

plt.plot(YPHI[:,9]*100000,XPHI[:,9]*100000,'ok',label ="zone radius 1")
plt.plot(YPHI[:,8]*100000,XPHI[:,8]*100000,'xk',label="zone radius 0.9")
plt.plot(YPHI[:,7]*100000,XPHI[:,7]*100000,'+k',label="zone radius 0.8")
plt.plot(YPHI[:,6]*100000,XPHI[:,6]*100000,'^k',label="zone radius 0.7")
plt.plot(YPHI[:,5]*100000,XPHI[:,5]*100000,'sk',label="zone radius 0.6")
plt.plot(YPHI[:,0]*100000,XPHI[:,0]*100000,'sk',label="zone radius 0.1")
plt.xlabel("Y (image space) x100k",fontsize=10)
plt.ylabel("X (image space) x100k",fontsize=10)
plt.xticks(fontsize = 7)
plt.yticks(fontsize = 7)
#plt.legend(loc="upper left")
plt.axis('scaled')

plt.show()
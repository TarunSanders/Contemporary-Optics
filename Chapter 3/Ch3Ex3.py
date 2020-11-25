# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:21:55 2020

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

#definition of rays
cx1top = 0.01
cy1top = 0
cz1top = math.sqrt(1-cx1top**2-cy1top**2)
#cz1top = math.sqrt(1-cx1top**2-cy1top**2)

x1top = np.arange(1,0,-0.1)
y1top = 0
FP = np.zeros(x1top.shape)
#z1top = r1 - math.sqrt(r1**2 - x1top**2)

cx1bot = cx1top
cy1bot = 0
cz1bot = cz1top
# cz1bot = math.sqrt(1-cx1top**2-cy1top**2)

# x1bot = -x1top
y1bot = y1top

# z1bot = r1 - math.sqrt(r1**2 - x1bot**2)

Phi = math.pi/180*np.arange(0,180,15)

XPHI = []
YPHI = []

for j, x1T in enumerate(x1top):
#emerging top ray parameters
  
  
    z1top = r1 - math.sqrt(r1**2 - x1T**2)
   
    x1B = -x1T
    z1bot = z1top
    x2dashtop,y2dashtop,z2dashtop, cx2dashtop, cy2dashtop,cz2dashtop = rayLens(x1T,y1top,z1top,cx1top,cy1top,cz1top,n1,n1dash,t,r1,r2)
     
    x2dashbot,y2dashbot,z2dashbot, cx2dashbot, cy2dashbot,cz2dashbot = rayLens(x1B,y1bot,z1bot,cx1bot,cy1bot,cz1bot,n1,n1dash,t,r1,r2)

    FP[j] = findIntersection(x2dashtop, x2dashbot, z2dashtop, z2dashbot, cx2dashtop, cx2dashbot, cz2dashtop, cz2dashbot)
    


       #xFP,yFP,zFP = rayFPplane(FP, x2dashtop, y2dashtop, z2dashtop, cx2dashtop, cy2dashtop, cz2dashtop, n1)
#for ii in range(0,x1top.shape[0]-1):
   
for ii in range(0,x1top.shape[0]-1):
    x3ar = []
    y3ar = []
    for jj in range(0,Phi.shape[0]-1):
        x1 = x1top[ii]*math.cos(Phi[jj])
        y1 = x1top[ii]*math.sin(Phi[jj])
        z1 = r1 - math.sqrt(r1**2 - x1top[ii]**2)
        
        x2dash,y2dash,z2dash, cx2dash, cy2dash,cz2dash = rayLens(x1,y1,z1,cx1top,cy1top,cz1top,n1,n1dash,t,r1,r2)
        x3,y3,z3 = rayFPplane(FP[ii],x2dash,y2dash,z2dash, cx2dash, cy2dash,cz2dash,n1)
        x3ar.append(x3)
        y3ar.append(y3)
    XPHI.append(x3ar)
    YPHI.append(y3ar)
# XPHI.reshape(-1,x1top.shape[0])
# YPHI.reshape(-1,x1top.shape[0])

plt.plot(YPHI[0][:],XPHI[0][:],'ok',label = "zone radius 1")
plt.plot(YPHI[1][:],XPHI[1][:],'xk',label="zone radius 0.9")
plt.plot(YPHI[2][:],XPHI[2][:],'+k',label="zone radius 0.8")
plt.plot(YPHI[3][:],XPHI[3][:],'^k',label="zone radius 0.7")
plt.plot(YPHI[4][:],XPHI[4][:],'sk',label="zone radius 0.6")

#plt.xaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel("Y (image space)",fontsize=20)
plt.ylabel("X (image space)",fontsize=20)
plt.xticks(fontsize = 7)
plt.yticks(fontsize = 7)
plt.legend(loc="upper left")
plt.axis('scaled')
plt.xlim(-0.00055,0.00055)
plt.ylim(0.10225,0.1042)

plt.show()
#print(FP)
Phi = math.pi/180*np.arange(180,361,15)
XPHIbot = []
YPHIbot = []
for ii in range(0,x1top.shape[0]-1):
    x3ar = []
    y3ar = []
    for jj in range(0,Phi.shape[0]-1):
        x1 = x1top[ii]*math.cos(Phi[jj])
        y1 = x1top[ii]*math.sin(Phi[jj])
        z1 = r1 - math.sqrt(r1**2 - x1top[ii]**2)
        
        x2dash,y2dash,z2dash, cx2dash, cy2dash,cz2dash = rayLens(x1,y1,z1,cx1top,cy1top,cz1top,n1,n1dash,t,r1,r2)
        x3,y3,z3 = rayFPplane(FP[ii],x2dash,y2dash,z2dash, cx2dash, cy2dash,cz2dash,n1)
        x3ar.append(x3)
        y3ar.append(y3)
    XPHIbot.append(x3ar)
    YPHIbot.append(y3ar)
plt.plot(YPHI[0][:],XPHI[0][:],'ok',label = "0 to 180 deg")
plt.plot(YPHIbot[0][:],XPHIbot[0][:],'or',label = "180 to 360 deg")
plt.xlabel("Y (image space)",fontsize=20)
plt.ylabel("X (image space)",fontsize=20)
plt.xticks(fontsize = 7)
plt.yticks(fontsize = 7)
plt.legend(loc="upper left")
plt.axis('scaled')
plt.xlim(-0.00055,0.00055)
plt.ylim(0.10225,0.1038)
#plt.xaxis.set_major_locator(plt.MaxNLocator(2))
plt.show()
#plt.xlim(-2,2)
#plt.ylim(0,-1.2)

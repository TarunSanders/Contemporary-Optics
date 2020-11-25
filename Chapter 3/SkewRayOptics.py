# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 14:00:56 2020

@author: LabUser
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

#inputs: coordinates of object, position of object from vertex of lens, directional cosines of ray with axes, radius of curvature of surface, 
def cosTheta(x,y,z,lp, radius,cx,cy,cz):
    Sigma = x*cx+y*cy+z*cz
    L = z - lp #distance of first origin from vertex
    c = L + radius # distance of origin to centre of curvature
    Curv = 1/radius  #curvature
    cosTheta = math.sqrt((Curv**2)*((cz*c-Sigma)**2)-(Curv**2)*(x**2+y**2+(z-L)**2)+2*Curv*(z-L))
    return cosTheta

def rayDisplacement(x,y,z,lp,radius,cx,cy,cz,cosTheta):
    Sigma = x*cx+y*cy+z*cz
    Curv = 1/radius
    L = z - lp
    c = L + radius
    
    T1 = (Curv*(cz*c-Sigma)-cosTheta)/Curv
    
    return T1

def skewPower(n1,n1Dash,radius,cosTheta):
    Theta = math.acos(cosTheta) #incident angle
    ThetaDash = math.asin(n1*math.sin(Theta)/n1Dash)
    cosThetaDash = math.cos(ThetaDash)
    K1 = (n1Dash*cosThetaDash-n1*cosTheta)/radius
    return K1

# for z equation : increases z to z1+L as the z-cordinate at point of incidence is taken with respect to the  vertex V1
def TranslationMatrix(n1,T1):
    T = np.array([[1,0],[T1/n1,1]])
    return T

# for z equation specifically the n1Dash*czDash - K1/Curv = n1*cz - K1*z1
def RefractionMatrix(K1):
    R = np.array([[1,-K1],[0,1]])
    return R
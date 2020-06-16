#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 22:13:16 2020

@author: youssef
"""

import numpy as np
import scipy.interpolate as sp
import matplotlib.pyplot as plt

Altitude = np.arange(0,13000,1000)
Mach     = np.array([0,0.25,0.5,0.78])
Values   = np.array([[96722.5216	,   67428.05362,47797.17295	,34848.43582],
                     [92730.66919,  66953.61247, 49363.17223,34875.76088],
                     [83327.33083,	60272.56754,	44580.51898,	31811.88102],
                     [75528.07162,	54984.87745,	41105.11362,	30166.20097],
                     [68192.30809,	49935.16789,	37672.74605,	28178.95313],
                     [61325.46673,	45143.27947,	34329.63267,	26076.81658],
                     [54934.35083,	40632.33221,	31115.11061,	23940.46769],
                     [52728.40771,	38788.01903,	29339.68686,	21994.39766],
                     [47022.7086,	34789.833,	26546.78776,	20260.27162],
                     [41496.00717,	30791.69771,	23600.46648,	18162.6575],
                     [36575.19359,	27245.49901,	21000.56856,	16329.92886],
                     [33957.59999,	25411.71708,	19713.45868,	15500.42165],
                     [27471.22224,	20554.11247,	15941.42122,	12529.85862]])

Thrust = sp.interp2d(Mach,Altitude,Values,kind='cubic')

mach     = np.arange(0,0.71,0.1)
thrust   = Thrust(mach,Altitude)


# TSFC function with interpolation
Altitude = np.arange(0,16000,4000)
Mach     = np.arange(0,1,0.5)
Values   = (10**-6)*np.array([[5.1996, 10.07987],
                     [5.53191581,	9.747087486],
                     [6.117945716,	10.29400869],
                     [6.655300979,	10.95611142]])


TSFCf    = sp.interp2d(Mach,Altitude,Values,kind='linear')
TSFC = TSFCf(Mach,Altitude)







plt.plot(Altitude,TSFC[:,0])
plt.plot(Altitude,TSFC[:,1])
plt.scatter(Altitude,Values[:,0])
plt.scatter(Altitude,Values[:,1])
plt.xlabel("Altitude")
plt.ylabel("TSFC")
plt.grid()
plt.show()
# -*- coding: utf-8 -*-
"""

"""


import numpy as np
from matplotlib import pyplot as plt



N = 100
delta = 3.83e-1
delta_b = 1.24
V_inf = 230.13 

y = np.linspace(0,delta,N)
y_b = np.linspace(0,delta_b,N)
vel = V_inf*((y/delta)**(1/7))
vel_b = V_inf*((y_b/delta_b)**(1/7))

plt.figure()
ax = plt.axes()
# plt.vlines(230.13,0,delta,color="r",label="V Free Stream")
plt.plot(vel,y)
plt.xlabel("Velocity [m/s]")
plt.ylabel("y [m]")
plt.legend()
plt.grid()
plt.title("Top Surface (x = 41.6 m)")



plt.figure()
plt.plot(vel_b,y_b)
# plt.vlines(230.13,0,delta_b,color="r",label="V Free Stream")
plt.xlabel("Velocity [m/s]")
plt.ylabel("y [m]")
plt.legend()
plt.grid()
plt.title("Bottom Surface (x = 41.6 m)")
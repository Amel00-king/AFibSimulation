# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:51:06 2023

@author: Edu
"""

# FitzHugh-Nagumo (FN) equations solver 
# Libraries:
import numpy as np
import pylab as plt

# Size of mesh:
S = 100
# Create grid of 100 by 100:
X,Y = np.meshgrid(np.arange(S), np.arange(S))

sze = (S,S)
Ugrid = np.zeros(sze)
Ngrid = np.zeros(sze)

# Set initial conditions:
for m in range(0,3):
    for n in range(0,S):
        Ugrid[n,m] = 1.0

# Numerical values:
D = 0
eps = 0.02
g = 5.0
k = 8.0
u_th = 0.1
dt = 0.1
dx = 0.2
dy = 0.2
D = 0.1



# Central differences function for 2nd pd:
def CD2(val1, val2, val3, dw):
    deriv2nd  = (val1 - 2*val2 + val3) / (dw*dw)
    return deriv2nd


# Forward Euler method and Central Differences for solving:
for v in range(0,1000):
    Ugrid2 = np.zeros(sze)
    Ngrid2 = np.zeros(sze)
    for i in range(0,S-1):
        for j in range(0,S-1):
            Ugrid2[i,j] = Ugrid[i,j] + dt*(D*(CD2(Ugrid[i+1,j],Ugrid[i,j],Ugrid[i-1,j], dx) +
                                              CD2(Ugrid[i,j+1],Ugrid[i,j],Ugrid[i,j-1], dy)) - (k*Ugrid[i,j]*(Ugrid[i,j]-u_th)*(Ugrid[i,j]-1)))
                                              
            Ngrid2[i,j] = Ngrid[i,j] + dt*(eps * (g*Ugrid[i,j] - Ngrid[i,j]))
            
            # Boundry conditions are derivative are same at border:
            for b in range(0,S):
                Ugrid2[b,0] = Ugrid2[b,1]
                Ugrid2[b,S-1] = Ugrid2[b, S-2]
                Ugrid2[0,b] = Ugrid2[1,b]
                Ugrid2[S-1,b] = Ugrid2[S-2,b]
    
    # Reset old grids:
    for i in range(0,S):
        for j in range(0,S):
                Ugrid[i,j] = Ugrid2[i,j]
                Ngrid[i,j] = Ngrid2[i,j]
    
    if (v % 10) == 0:
        plt.pcolormesh(X,Y,Ugrid)
        plt.show()
    
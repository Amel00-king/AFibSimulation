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
dt = 0.05
dx = 0.2
dy = 0.2
D = 0.1

# Central differences function for 2nd pd:
def CD2(val1, val2, val3, dw):
    deriv2nd  = (val1 - 2*val2 + val3) / (dw*dw)
    return deriv2nd


# Forward Euler method and Central Differences for solving:
only_once = True
for itern in range(0,2000):
    print('Iteration, ', itern)
    Ugrid2 = np.zeros(sze)
    Ngrid2 = np.zeros(sze)
    for i in range(0,S-1):            
        for j in range(0,S-1):
            Ugrid2[i,j] = Ugrid[i,j] + dt*(D*(CD2(Ugrid[i+1,j],Ugrid[i,j],Ugrid[i-1,j], dx) +
                                                  CD2(Ugrid[i,j+1],Ugrid[i,j],Ugrid[i,j-1], dy)) - (k*Ugrid[i,j]*(Ugrid[i,j]-u_th)*(Ugrid[i,j]-1)) - Ngrid[i,j])
                                                  
            Ngrid2[i,j] = Ngrid[i,j] + dt*(eps * (g*Ugrid[i,j] - Ngrid[i,j]))
            
    # Make spiral wave set all vals to 0:
    if ((Ugrid2[int(S/2),int(S/2)] > 0.2) & (only_once)):
        for rep in range(int(S/2), S):
            for rep2 in range(0,S):
                Ugrid2[rep,rep2] = 0
                #Ngrid2[rep,rep2] = 0
        only_once = False
        print('Completed replacement')
        
    # Make fibrosis patch for extra signal:
    if (itern == 1000):
        for patch in range(15, 20):
            for patch2 in range(15,20):
                Ugrid2[patch,patch2] = 1
                Ngrid2[patch,patch2] = 0
        
                
    # Boundry conditions are derivative is 0 at border:
    for b in range(0,S):
        Ugrid2[b,0] = Ugrid2[b,1]
        Ugrid2[b,S-1] = Ugrid2[b, S-2]
        Ugrid2[0,b] = Ugrid2[1,b]
        Ugrid2[S-1,b] = Ugrid2[S-2,b]
    
    # Reset old grids:
    Ugrid = Ugrid2
    Ngrid = Ngrid2
    
    if (itern % 10) == 0:
        plt.pcolormesh(X,Y,Ugrid)
        plt.show()
        plt.title(itern)
    
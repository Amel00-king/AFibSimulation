# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:41:30 2023

@author: olegadmin
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
eps = 0.02
g = 5.0
k = 8.0
u_th = 0.1
dt = 0.05
dx = 0.2
dy = 0.2
D = 0.1

# Values for the numer of subdivisions:
# If both are 2 then only one spiral in center will form.
n_divs = 2
m_divs = 2

# Value for the electrode placements:
# Assume all electrode placements in a line across atrium **DIDNT WORK
# Place electrodes according to how it looks in pictures, taking MESH SIZE as 100
R = [150,115]
L = [-70,115]
F = [-15,-15]
L1 = [-15,-15]
L2 = [80,-15]
L3 = [90,-20]
L4 = [100,-25]
L5 = [110,-30]
L6 = [110,-30]
xe = 5
ye = 5
ze = 40
ze6 = 40

# Visualise electrode placement:
plt.plot([R[0],L[0],F[0],L1[0],L2[0],L2[0],L3[0],L4[0],L5[0],L6[0]], [R[1],L[1],F[1],L1[1],L2[1],L2[1],L3[1],L4[1],L5[1],L6[1]], 'ro')
plt.Rectangle((0, 0), 100, 100)
plt.xlim([-100,200])
plt.ylim([-50,150])
plt.title('Electrode placement shown by red dot:')
plt.show()

# Central differences function for 2nd pd:
def CD2(val1, val2, val3, dw):
    deriv2nd  = (val1 - 2*val2 + val3) / (dw*dw)
    return deriv2nd                                        

# Place fibrosis patch:                                    
def PlaceFib(Grid1, Grid2, min_y, max_y, min_x, max_x):
    for patch in range(min_y, max_y):
        for patch2 in range(min_x, max_x):
            Grid1[patch,patch2] = 1
            Grid2[patch,patch2] = 0
    return Grid1, Grid2

# Get focal points and spiral creation points:
# Returns 2d array of points where focal and spiral/cut-off points are
def GetArray(Grid, n, m):
    array = np.zeros((int((n-1)*(m-1)), 2), dtype = int)
    multx = 1
    sze = Grid.shape
    x_step = int(sze[0] / n)
    y_step = int(sze[1] / m)
    iteration = 0
    for i in range(0,n-1):
        multy = 1
        for j in range(0,m-1):
            array[iteration,0] = x_step*multx
            array[iteration,1] = y_step*multy
            multy += 1
            iteration += 1
        multx += 1
    return array

# Spiral wave generation:
def MakeSpiral(Grid, x_val, S):
    for rep in range(x_val, S):
        for rep2 in range(0,S):
            Grid[rep,rep2] = 0
    only_once = False
    print('Completed replacement')
    return Grid, only_once
        

# Get array for spira/focal point coordinates:
Coord_array = GetArray(Ugrid, n_divs, m_divs)

# Forward Euler method and Central Differences for solving:
tot_iterations = 20

# Initialise set of all P-wave traces for each electrode lead:
All_Pwaves_Lii = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_Liii = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_WCT = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_V1 = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_V2 = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_V3 = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_V4 = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_V5 = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_V6 = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)
All_Pwaves_aVF = np.zeros((((n_divs-1)*(m_divs-1)), tot_iterations), dtype=float)

# Loop through for all possible values for the spiral location intiation:
for spiral_pos in range(0,((n_divs-1)*(m_divs-1))):
    print('Position of spiral', (spiral_pos+1), ' is: ', Coord_array[spiral_pos])
    
    # Create grid of 100 by 100:
    sze = (S,S)
    Ugrid = np.zeros(sze)
    Ngrid = np.zeros(sze)
    
    # Set initial conditions:
    for m in range(0,3):
        for n in range(0,S):
            Ugrid[n,m] = 1.0
    
    # Condition for spiral wave set only once:
    only_once = True
    
    # Initialise array for ECG values:
    Int_ECG_R = np.zeros(tot_iterations)
    Int_ECG_L = np.zeros(tot_iterations)
    Int_ECG_F = np.zeros(tot_iterations)
    Int_ECG_L1 = np.zeros(tot_iterations)
    Int_ECG_L2 = np.zeros(tot_iterations)
    Int_ECG_L3 = np.zeros(tot_iterations)
    Int_ECG_L4 = np.zeros(tot_iterations)
    Int_ECG_L5 = np.zeros(tot_iterations)
    Int_ECG_L6 = np.zeros(tot_iterations)
    

        
    for itern in range(0,tot_iterations):
        #print('Iteration, ', itern)
        # Initialise update grids:
        Ugrid2 = np.zeros(sze)
        Ngrid2 = np.zeros(sze)
        sum_R = 0
        sum_L = 0
        sum_F = 0
        sum_L1 = 0
        sum_L2 = 0
        sum_L3 = 0
        sum_L4 = 0
        sum_L5 = 0
        sum_L6 = 0
        
        for i in range(0,S-1):            
            for j in range(0,S-1):
                Laplacian_ij = D*(CD2(Ugrid[i+1,j],Ugrid[i,j],Ugrid[i-1,j], dx) + CD2(Ugrid[i,j+1],Ugrid[i,j],Ugrid[i,j-1], dy))
                Ugrid2[i,j] = Ugrid[i,j] + dt*(Laplacian_ij - (k*Ugrid[i,j]*(Ugrid[i,j]-u_th)*(Ugrid[i,j]-1)) - Ngrid[i,j])
                                                      
                Ngrid2[i,j] = Ngrid[i,j] + dt*(eps * (g*Ugrid[i,j] - Ngrid[i,j]))#
                sum_R += Laplacian_ij / pow((pow((i-R[0]),2) + pow((j-R[1]), 2) + pow(ze,2)), 0.5)
                sum_L += Laplacian_ij / pow((pow((i-L[0]),2) + pow((j-L[1]), 2) + pow(ze,2)), 0.5)
                sum_F += Laplacian_ij / pow((pow((i-F[0]),2) + pow((j-F[1]), 2) + pow(ze,2)), 0.5)
                sum_L1 += Laplacian_ij / pow((pow((i-L1[0]),2) + pow((j-L1[1]), 2) + pow(ze,2)), 0.5)
                sum_L2 += Laplacian_ij / pow((pow((i-L2[0]),2) + pow((j-L2[1]), 2) + pow(ze,2)), 0.5)
                sum_L3 += Laplacian_ij / pow((pow((i-L3[0]),2) + pow((j-L3[1]), 2) + pow(ze,2)), 0.5)
                sum_L4 += Laplacian_ij / pow((pow((i-L4[0]),2) + pow((j-L4[1]), 2) + pow(ze,2)), 0.5)
                sum_L5 += Laplacian_ij / pow((pow((i-L5[0]),2) + pow((j-L5[1]), 2) + pow(ze,2)), 0.5)
                sum_L6 += Laplacian_ij / pow((pow((i-L6[0]),2) + pow((j-L6[1]), 2) + pow(ze,2)), 0.5)
                
        # Make spiral wave set all vals to 0:
        # Spiral location is the first value in the Coord_grid array:
        location = Coord_array[spiral_pos]
        if ((Ugrid2[location[0],location[1]] > 0.2) & (only_once)):
            Ugrid2, only_once = MakeSpiral(Ugrid2, location[0], S)
         
        # NOT NECESSARY CURRENTLY:
        # Make fibrosis patch for extra signal:
        #if (itern == 800):
        #    Ugrid2, Ngrid2 = PlaceFib(Ugrid2, Ngrid2, 15, 25, 15, 25)
        
    
        # Boundry conditions are derivative is 0 at border:
        for b in range(0,S):
            Ugrid2[b,0] = Ugrid2[b,1]
            Ugrid2[b,S-1] = Ugrid2[b, S-2]
            Ugrid2[0,b] = Ugrid2[1,b]
            Ugrid2[S-1,b] = Ugrid2[S-2,b]
            
        # Calculate integral for Int_ECG for each electrode placement:
        # Sum of all laplacians:
        Int_ECG_R[itern] = round(sum_R,4)
        Int_ECG_L[itern] = round(sum_L,4)
        Int_ECG_F[itern] = round(sum_F,4)
        Int_ECG_L1[itern] = round(sum_L1,4)
        Int_ECG_L2[itern] = round(sum_L2,4)
        Int_ECG_L3[itern] = round(sum_L3,4)
        Int_ECG_L4[itern] = round(sum_L4,4)
        Int_ECG_L5[itern] = round(sum_L5,4)
        Int_ECG_L6[itern] = round(sum_L6,4)
            
        # Reset old grids:
        Ugrid = Ugrid2
        Ngrid = Ngrid2
        
        # Plot to show updating mesh:
        if (itern % 10) == 0:
            plt.pcolormesh(X,Y,Ugrid)
            plt.title(itern)
            plt.show()
        print('Completed iteration: ', itern)
    
    # Add each lead ECG wave to corresponding ECG array set:
    All_Pwaves_Lii[spiral_pos] = Int_ECG_F - Int_ECG_R
    All_Pwaves_Liii[spiral_pos] = Int_ECG_F - Int_ECG_L
    WCT = (Int_ECG_F + Int_ECG_R + Int_ECG_L) / 3
    All_Pwaves_WCT[spiral_pos] = WCT
    All_Pwaves_V1[spiral_pos] = Int_ECG_L1 - WCT
    All_Pwaves_V2[spiral_pos] = Int_ECG_L2 - WCT
    All_Pwaves_V3[spiral_pos] = Int_ECG_L3 - WCT
    All_Pwaves_V4[spiral_pos] = Int_ECG_L4 - WCT
    All_Pwaves_V5[spiral_pos] = Int_ECG_L5 - WCT
    All_Pwaves_V6[spiral_pos] = Int_ECG_L6 - WCT
    All_Pwaves_aVF[spiral_pos] = Int_ECG_F - ((Int_ECG_R+Int_ECG_L)/2)


# Plot all ECG traces from all 9 electrodes **ASSUMING ONLY 1 SPIRAL**
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_Lii), label='Lii')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_Liii), label='Liii')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_WCT), label='WCT')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_V1), label='V1')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_V2), label='V2')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_V3), label='V3')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_V4), label='V4')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_V5), label='V5')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_V6), label='V6')
plt.plot(np.arange(tot_iterations), np.transpose(All_Pwaves_aVF), label='aVF')
plt.legend(loc=0, prop={'size': 6})
plt.title('All 9 electrode ECGs: ')
plt.show()

# Loop through all spiral locations and set each ECG value for columnisation:
for i in range(0,(n_divs-1)*(m_divs-1)):
    All_text_array = np.around(np.stack((All_Pwaves_Lii[i], All_Pwaves_Liii[i], All_Pwaves_WCT[i], All_Pwaves_V1[i], All_Pwaves_V2[i], All_Pwaves_V3[i], All_Pwaves_V4[i], All_Pwaves_V5[i], All_Pwaves_V6[i], All_Pwaves_aVF[i]), axis=-1), decimals=4)
    filename = 'ECG9trace_' + str(i) + '_iterations_' + str(tot_iterations) + '_position_' + str(Coord_array[i][0]) + 'comma' + str(Coord_array[i][1]) + '_LiiLiiiWTCV1V2V3V4V5V6aVF.txt'
    with open(filename, 'w') as f:
        np.savetxt(filename, All_text_array, fmt='%1.4f', delimiter="\t")

# Display values/settings for screenshots:
print('Size of mesh = ', S, ' by ', S)
print('Number of spiral points = ', ((n_divs-1)*(m_divs-1)))
print('Spiral point coordinates: ', Coord_array)
print('Total iterations = ', tot_iterations)

# Method 1:
# 5 by 5 points each with ECG trace
# focal and spiral points

# FitzHugh-Nagumo (FN) equations solver for SPIRALS
# Libraries:
import numpy as np
import pylab as plt

# PARAMETERS:
# Size of mesh:
S = 100
# Create grid of 100 by 100:
X,Y = np.meshgrid(np.arange(S), np.arange(S))

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
n_dim = 10
m_dim = 10

# Define centers:
Center1 = [30,75]
Center2 = [75,30]

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

sze = (S,S)
Ugrid = np.zeros(sze)
Ngrid = np.zeros(sze)
#-----------------------------------

# FUNCTIONS:
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
def GetArray(center, n_dim, m_dim):
    array = np.zeros((int(n_dim*m_dim), 2), dtype = int)
    iteration = 0
    for i in range(-int(n_dim/2), int(n_dim/2)):
        for j in range(-int(m_dim/2), int(m_dim/2)):
            array[iteration,0] = center[0] + i
            array[iteration,1] = center[1] + j
            iteration +=1
    return array

# Spiral wave generation:
def MakeSpiral(Grid, x_val, S):
    for rep in range(x_val, S):
        for rep2 in range(0,S):
            Grid[rep,rep2] = 0
    only_once = False
    print('Completed replacement')
    return Grid, only_once
        
# TIC/TOC function - not original code
def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set"  )     

tic()
# Get array for spira/focal point coordinates:
Coord_array1 = GetArray(Center1, n_dim, m_dim)
Coord_array2 = GetArray(Center2, n_dim, m_dim)

All_Centers = [Center1, Center2]
All_pos_set = [Coord_array1, Coord_array2]

#----------------------------------
# FOR SPIRALS:
# Set initial conditions:
for m in range(0,3):
    for n in range(0,S):
        Ugrid[n,m] = 1.0

# Forward Euler method and Central Differences for solving:
tot_iterations = 6000
# When ECGs are stable so begin recording:
start_recording = 2500

eff_iterns = tot_iterations-start_recording

# Initialise array for ECG values:
Int_ECG_R = np.zeros(eff_iterns)
Int_ECG_L = np.zeros(eff_iterns)
Int_ECG_F = np.zeros(eff_iterns)
Int_ECG_L1 = np.zeros(eff_iterns)
Int_ECG_L2 = np.zeros(eff_iterns)
Int_ECG_L3 = np.zeros(eff_iterns)
Int_ECG_L4 = np.zeros(eff_iterns)
Int_ECG_L5 = np.zeros(eff_iterns)
Int_ECG_L6 = np.zeros(eff_iterns)

for itr in range(len(All_Centers)):
    Center_iter = All_Centers[itr]
    Coord_array = All_pos_set[itr]
    print('SPIRALS:')
    print('Center is ', Center_iter)
    print('All initiation positions: ', Coord_array)
    # Loop through for all possible values for the SPIRAL LOCATIONS:
    for spiral_pos in range(0,(n_dim*m_dim)):
        print('Position of SPIRAL ', (spiral_pos+1), ' is: ', Coord_array[spiral_pos])
        
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
                                                          
                    Ngrid2[i,j] = Ngrid[i,j] + dt*(eps * (g*Ugrid[i,j] - Ngrid[i,j]))
                    if (itern > start_recording):
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
            
        
            # Boundry conditions are derivative is 0 at border:
            for b in range(0,S):
                Ugrid2[b,0] = Ugrid2[b,1]
                Ugrid2[b,S-1] = Ugrid2[b, S-2]
                Ugrid2[0,b] = Ugrid2[1,b]
                Ugrid2[S-1,b] = Ugrid2[S-2,b]
                
            # Calculate integral for Int_ECG for each electrode placement:
            # Sum of all laplacians:
            if (itern > start_recording):
                Int_ECG_R[itern-start_recording] = round(sum_R,4)
                Int_ECG_L[itern-start_recording] = round(sum_L,4)
                Int_ECG_F[itern-start_recording] = round(sum_F,4)
                Int_ECG_L1[itern-start_recording] = round(sum_L1,4)
                Int_ECG_L2[itern-start_recording] = round(sum_L2,4)
                Int_ECG_L3[itern-start_recording] = round(sum_L3,4)
                Int_ECG_L4[itern-start_recording] = round(sum_L4,4)
                Int_ECG_L5[itern-start_recording] = round(sum_L5,4)
                Int_ECG_L6[itern-start_recording] = round(sum_L6,4)
                
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
        Lii = Int_ECG_F - Int_ECG_R
        Liii = Int_ECG_F - Int_ECG_L
        WCT = (Int_ECG_F + Int_ECG_R + Int_ECG_L) / 3
        V1 = Int_ECG_L1 - WCT
        V2 = Int_ECG_L2 - WCT
        V3 = Int_ECG_L3 - WCT
        V4 = Int_ECG_L4 - WCT
        V5 = Int_ECG_L5 - WCT
        V6 = Int_ECG_L6 - WCT
        aVF = Int_ECG_F - ((Int_ECG_R+Int_ECG_L)/2)
    
        # For plots see situation specific solvers/code
        
        All_text_array = np.around(np.stack((Lii, Liii, WCT, V1, V2, V3, V4, V5, V6, aVF), axis=-1), decimals=4)
        filename = 'ECGtrace_SP_' + str(spiral_pos) +'_Center_' + str(Center_iter[0]) + '_' + str(Center_iter[1]) + '_iterns_' + str(tot_iterations) + '_position_' + str(Coord_array[spiral_pos][0]) + 'comma' + str(Coord_array[spiral_pos][1]) + '_LiiLiiiWTCV1V2V3V4V5V6aVF.txt'
        with open(filename, 'w') as f:
            np.savetxt(filename, All_text_array, fmt='%1.4f', delimiter="\t")

# Display values/settings for screenshots:
print('Size of mesh = ', S, ' by ', S)
print('Number of spiral points = ', (n_dim*m_dim))
print('Spiral point coordinates: ', Coord_array[spiral_pos])
print('Total iterations = ', tot_iterations)

# FOR FOCAL POINTS:
# Size of focal point (6x6)
Focal_size = 6

# Forward Euler method and Central Differences for solving:
tot_iterations = 4000
# When ECGs are stable so begin recording:
start_recording = 500

eff_iterns = tot_iterations-start_recording

# Initialise array for ECG values:
Int_ECG_R = np.zeros(eff_iterns)
Int_ECG_L = np.zeros(eff_iterns)
Int_ECG_F = np.zeros(eff_iterns)
Int_ECG_L1 = np.zeros(eff_iterns)
Int_ECG_L2 = np.zeros(eff_iterns)
Int_ECG_L3 = np.zeros(eff_iterns)
Int_ECG_L4 = np.zeros(eff_iterns)
Int_ECG_L5 = np.zeros(eff_iterns)
Int_ECG_L6 = np.zeros(eff_iterns)

for itr in range(len(All_Centers)):
    Center_iter = All_Centers[itr]
    Coord_array = All_pos_set[itr]
    print('FOCAL POINTS:')
    print('Center is ', Center_iter)
    print('All initiation positions: ', Coord_array)
    # Loop through for all possible values for the FOCAL POINT LOCATIONS:
    for spiral_pos in range(0,(n_dim*m_dim)):
        print('Position of FOCAL POINT ', (spiral_pos+1), ' is: ', Coord_array[spiral_pos])
        
        # Create grid of 100 by 100:
        sze = (S,S)
        Ugrid = np.zeros(sze)
        Ngrid = np.zeros(sze)
        
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
            
            # Place patches every 1000 iterations:
            if (itern == 0) | (itern%970 == 0):
                Ugrid2, Ngrid2 = PlaceFib(Ugrid, Ngrid, (Coord_array[spiral_pos][0]-int(Focal_size / 2)), (Coord_array[spiral_pos][0]+int(Focal_size / 2)), (Coord_array[spiral_pos][1]-int(Focal_size / 2)),
                                          (Coord_array[spiral_pos][1]+int(Focal_size / 2)))
    
            for i in range(0,S-1):            
                for j in range(0,S-1):
                    Laplacian_ij = D*(CD2(Ugrid[i+1,j],Ugrid[i,j],Ugrid[i-1,j], dx) + CD2(Ugrid[i,j+1],Ugrid[i,j],Ugrid[i,j-1], dy))
                    Ugrid2[i,j] = Ugrid[i,j] + dt*(Laplacian_ij - (k*Ugrid[i,j]*(Ugrid[i,j]-u_th)*(Ugrid[i,j]-1)) - Ngrid[i,j])
                                                          
                    Ngrid2[i,j] = Ngrid[i,j] + dt*(eps * (g*Ugrid[i,j] - Ngrid[i,j]))
                    if (itern > start_recording):
                        sum_R += Laplacian_ij / pow((pow((i-R[0]),2) + pow((j-R[1]), 2) + pow(ze,2)), 0.5)
                        sum_L += Laplacian_ij / pow((pow((i-L[0]),2) + pow((j-L[1]), 2) + pow(ze,2)), 0.5)
                        sum_F += Laplacian_ij / pow((pow((i-F[0]),2) + pow((j-F[1]), 2) + pow(ze,2)), 0.5)
                        sum_L1 += Laplacian_ij / pow((pow((i-L1[0]),2) + pow((j-L1[1]), 2) + pow(ze,2)), 0.5)
                        sum_L2 += Laplacian_ij / pow((pow((i-L2[0]),2) + pow((j-L2[1]), 2) + pow(ze,2)), 0.5)
                        sum_L3 += Laplacian_ij / pow((pow((i-L3[0]),2) + pow((j-L3[1]), 2) + pow(ze,2)), 0.5)
                        sum_L4 += Laplacian_ij / pow((pow((i-L4[0]),2) + pow((j-L4[1]), 2) + pow(ze,2)), 0.5)
                        sum_L5 += Laplacian_ij / pow((pow((i-L5[0]),2) + pow((j-L5[1]), 2) + pow(ze,2)), 0.5)
                        sum_L6 += Laplacian_ij / pow((pow((i-L6[0]),2) + pow((j-L6[1]), 2) + pow(ze,2)), 0.5)
        
            # Boundry conditions are derivative is 0 at border:
            for b in range(0,S):
                Ugrid2[b,0] = Ugrid2[b,1]
                Ugrid2[b,S-1] = Ugrid2[b, S-2]
                Ugrid2[0,b] = Ugrid2[1,b]
                Ugrid2[S-1,b] = Ugrid2[S-2,b]
                
            # Calculate integral for Int_ECG for each electrode placement:
            # Sum of all laplacians:
            if (itern > start_recording):
                Int_ECG_R[itern-start_recording] = round(sum_R,4)
                Int_ECG_L[itern-start_recording] = round(sum_L,4)
                Int_ECG_F[itern-start_recording] = round(sum_F,4)
                Int_ECG_L1[itern-start_recording] = round(sum_L1,4)
                Int_ECG_L2[itern-start_recording] = round(sum_L2,4)
                Int_ECG_L3[itern-start_recording] = round(sum_L3,4)
                Int_ECG_L4[itern-start_recording] = round(sum_L4,4)
                Int_ECG_L5[itern-start_recording] = round(sum_L5,4)
                Int_ECG_L6[itern-start_recording] = round(sum_L6,4)
                
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
        Lii = Int_ECG_F - Int_ECG_R
        Liii = Int_ECG_F - Int_ECG_L
        WCT = (Int_ECG_F + Int_ECG_R + Int_ECG_L) / 3
        V1 = Int_ECG_L1 - WCT
        V2 = Int_ECG_L2 - WCT
        V3 = Int_ECG_L3 - WCT
        V4 = Int_ECG_L4 - WCT
        V5 = Int_ECG_L5 - WCT
        V6 = Int_ECG_L6 - WCT
        aVF = Int_ECG_F - ((Int_ECG_R+Int_ECG_L)/2)
    
        # For plots see situation specific solvers/code
        
        All_text_array = np.around(np.stack((Lii, Liii, WCT, V1, V2, V3, V4, V5, V6, aVF), axis=-1), decimals=4)
        filename = 'ECGtrace_FP_' + str(spiral_pos) +'_Center_' + str(Center_iter[0]) + '_' + str(Center_iter[1]) + '_iterns_' + str(tot_iterations) + '_position_' + str(Coord_array[spiral_pos][0]) + 'comma' + str(Coord_array[spiral_pos][1]) + '_LiiLiiiWTCV1V2V3V4V5V6aVF.txt'
        with open(filename, 'w') as f:
            np.savetxt(filename, All_text_array, fmt='%1.4f', delimiter="\t")

    # Display values/settings for screenshots:
    print('Size of mesh = ', S, ' by ', S)
    print('Number of spiral points = ', (n_dim*m_dim))
    print('Focal point coordinates: ', Coord_array[spiral_pos])
    print('Total iterations = ', tot_iterations)

toc()
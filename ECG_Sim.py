'''
Running Fitz-Hugh Nagumo simulations on a 2D grid and the associated ECG signals 
'''
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd 



D = 0.5
epsilon = 0.02
g = 5.0
k = 8.0
uth = 0.1

#boundary conditions and steps
dx = 0.2
dy = 0.2
dt = 0.05 #used to be 0.1
L_x = 100
L_y = 100
frame_step = 100 #how many frames can pass before recording to a gif

#for rotor creation, reset up to the chosen height of the grid to 0 at chosen reset_times
#based on a 100x100 grid leaving space for the rotor to avoid washing out at the boundaries
reset_times = [950, 1425, 1900, 2376, 2850]
reset_heights = [30, 40, 50, 60, 70]

electrode_positions = {
    "R" : [150,115,40], #left shoulder
    "L" : [-70,115,40], #right shoulder
    "V1" : [-15,-15,40], 
    "V2" : [15,-15,40],
    "V3" : [45,-20,40],
    "V4" : [80,-25,40],
    "V5" : [110,-30,40],
    "V6" : [130,-35,40]
}

def plot_electrodes(electrodes):
    #plot the positions of the electrodes as well
    #extract x & y coordinates
    x_coords = [pos[0] for pos in electrodes.values()]
    y_coords = [pos[1] for pos in electrodes.values()]

    # electrode positionss
    plt.Rectangle((0, 0), 100, 100)
    plt.xlim([-100,200])
    plt.ylim([-50,150])
    plt.title('Electrode placement shown by red o:')
    plt.scatter(x_coords, y_coords, marker='o', color='red', label='Electrodes')
    plt.show()

    return 


''' only returns 100 values per electrode. It needs to calculate across x,y dimensionss not just one or other'''
def calc_distances(electrodes, L_x, L_y, dx, dy):
    distances = {}
    
    #calculate distance from each electrode to every element on the 2D 100x100 simulation grid
    for label, pos in electrodes.items():
        x, y, z = pos
        distance = np.sqrt((np.arange(L_x) - x)**2 + (np.arange(L_y) - y)**2 + (0 - z)**2)
        distances[label] = distance

    
    #returns a dictionary of keys with electrode names, and values as arrays of scalar float distances to grid
    return distances

def simulate(time_run, reset_time, reset_height, electrodes=None):
    u = np.zeros((L_x,L_y))
    n = np.zeros((L_x,L_y))
    u_history = [u]
    u_next = np.zeros((L_x,L_y))

    #initialize u for wave propogation
    u[:,:3] = 1
    
    # calculate the scalar distances for the grid to electrodes
    electrode_distances = calc_distances(electrode_positions,L_x,L_y, dx, dy)



    electrode_labels = list(electrodes.keys())
    electrode_data = {label: [] for label in electrode_labels}

    for t in range(time_run):

        if t == reset_time:
        print(f"creating rotor at: {t}")
        u[reset_height:, :] = 0
        u_history = [u]


        # Iterate through electrodes
        for label, distance in distances.items():
            x, y, z = electrodes[label]
            # Calculate integrand (replace the following lines with your actual calculations)
            sec_part_ux = np.random.rand(L_x, L_y)
            sec_part_uy = np.random.rand(L_x, L_y)
            integrand = (D * (sec_part_ux + sec_part_uy)) / distance

            # Trapezoidal rule integration
            integrated_signal = np.trapz(integrand)

            # Save the integrated signal for the current timestep
            electrode_data[label].append(integrated_signal)

        # Your other simulation logic here...

    # Convert electrode data to a DataFrame
    df = pd.DataFrame(electrode_data)

    # Plot electrode positions
    plot_electrode_positions(electrodes)

    return u_history, df



def simulate_all(electrode_pos, reset_times, reset_heights):

    for i in range(len(reset_times)):
        #It takes approx 4500-4600 time steps for a rotor wave to complete 1 cycle
        # therefore, total simulation time is 4700* 5 completed rotors, plus whatever time step the wave was initiated at
        sim_time = 4700 * 5 + reset_times(i)
        simulate(sim_time, reset_times(i), reset_heights(i), electrodes = electrode_pos)

    return  



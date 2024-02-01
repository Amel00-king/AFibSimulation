'''
Running Fitz-Hugh Nagumo simulations on a 2D grid and the associated ECG signals
'''
import numpy as np
import pandas as pd
#from google.colab import files, drive
import os
from tqdm import tqdm


D = 0.5
epsilon = 0.02
g = 5.0
k = 8.0
uth = 0.1

#boundary conditions and steps
dx = 0.2
dy = 0.2
dt = 0.01 #used to be 0.01
L_x = 100
L_y = 100
frame_step = 100 #how many frames can pass before recording to a gif

#for rotor creation, reset up to the chosen height of the grid to 0 at chosen reset_times
#based on a 100x100 grid leaving space for the rotor to avoid washing out at the boundaries
reset_times = [950, 1425, 1900, 2375, 2850]
reset_heights = [30, 40, 50, 60, 70]

electrode_positions = {
    "L" : [150,130,40], #left shoulder
    "R" : [-150,130,40], #right shoulder
    "V1" : [-15,15,40],
    "V2" : [15,15,40],
    "V3" : [45,0,40],
    "V4" : [80,-7,40],
    "V5" : [110,-15,40],
    "V6" : [130,-20,40],
    "LL" : [125,-150,40]
}



''' only returns 100 values per electrode --> attempted fix done'''
def calc_distances(electrodes, L_x, L_y, dx, dy, sim_time, frame_step):
    distance_arrays = []
    ecg_arrays = []

    for label, pos in electrodes.items():
        x, y, z = pos
        distance = np.zeros((L_y, L_x))

        #array of zeros to store ecg values every recorded frame_step
        ecg_values = np.zeros((int((sim_time/frame_step))+1))

        for i in range(L_x):
            for j in range(L_y):
                distance[j, i] = np.sqrt((i - x)**2 + (j - y)**2 + (0 - z)**2)

        distance_arrays.append(distance)
        ecg_arrays.append(ecg_values)

    #return counter for number of frames passed
    frame_count = 0

    return distance_arrays, ecg_arrays, frame_count


def simulate(time_run, reset_time, reset_height, electrodes=None):
    u = np.zeros((L_x,L_y))
    n = np.zeros((L_x,L_y))
    u_history = [u]
    u_next = np.zeros((L_x,L_y))

    #initialize u for wave propogation
    u[:,:3] = 1


    # Calculate the scalar distances for the grid to electrodes
    electrode_distances, ecg_signals, frame_count  = calc_distances(electrode_positions, L_x, L_y, dx, dy, time_run,frame_step)

    electrode_data = {label: [] for label in electrode_positions.keys()}

    for t in tqdm(range(time_run)):

        if t == reset_time:
            print(f"Creating rotor at: {t}")
            u[(reset_height):, :] = 0
            u_history = [u]

        for j in range(L_y -1):
            for i in range(L_x-1):
                
                current_element = u[i,j]
                n_elem = n[i,j]
                # Find central difference for each element
                CD_x = (u[i + 1, j] - 2 * current_element + u[i - 1, j]) / dx ** 2
                CD_y = (u[i, j + 1] - 2 * current_element + u[i, j - 1]) / dy ** 2

                
                # Update equations
                f = -k * current_element * (current_element - uth)  * (current_element - 1)
                n_elem = n_elem + dt * (epsilon * (g * n_elem - n_elem))
                u_next[i, j] = current_element + dt * (D * (CD_x + CD_y) + f - n_elem)

                # boundary conditions
                u_next[0,:] = u_next[1,:]
                u_next[L_x-1,:] = u_next[L_x-2,:]
                u_next[:,0] = u_next[:,1]
                u_next[:,L_x-1] = u_next[:,L_x-2]
                
                # Every x timesteps, record electrode signal
                if t % (frame_step )== 0:
                    # Iterate through electrodes
                    for label, electrode_dist, electrode_signal in zip(electrode_positions.keys(), electrode_distances, ecg_signals):

                        # Calculate integrand
                        integrand = (D * (CD_x + CD_y)) / electrode_dist[i, j]

                        # Add integrand of element to signal of entire grid in this timestep
                        electrode_signal[frame_count] += integrand



        #update u for the gif
        u = np.copy(u_next)

        if t % frame_step == 0:
            frame_count += 1
            u_history.append(np.copy(u))



    # Convert electrode data to a dictionary
    for label, electrode_signal in zip(electrode_positions.keys(), ecg_signals):
        electrode_data[label] = electrode_signal.tolist()

    # Plot electrode positions
    # plot_electrodes(electrodes)

    return u_history, electrode_data


def sim_all(reset_times, reset_heights, electrode_positions):

    current_dir = os.getcwd()

    for i in range(len(reset_times)):
        for j in range(len(reset_heights)):

            #loop through -2,-1,0,1,2 to get a 5x5 grid of the points
            for k in range(-2, 2):
                for l in range(-2, 2):

                    #sim_time = 300
                    sim_time = 4700 * 7 + reset_times[i]
                    # 38 timesteps per unit approx
                    reset = reset_times[i] + (k * 38)

                    height = reset_heights[j] + l

                    u_history,  electrode = simulate(sim_time, reset, height, electrode_positions)

                    df_electrode_data = pd.DataFrame(electrode)

                    # Save the DataFrame to a CSV file
                    if current_dir:
                        
                        csv_filename = os.path.join(current_dir, f'electrode_data_reset_{reset_times[i]}_height_{reset_heights[j]}_x_{k}_y_{l}.csv')
                        df_electrode_data.to_csv(csv_filename, index=False)
                        print(f"Electrode data saved to: {csv_filename}")

                    # update()
                    # #save the rotor gif to google drive too
                    # save_gif(u_history, electrode, gif_folder,reset_time=reset_times[i], reset_height=reset_heights[j])



sim_all(reset_times,reset_heights,electrode_positions)

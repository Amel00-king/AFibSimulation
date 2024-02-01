import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from numba import jit, njit, prange, set_num_threads

D = 0.5
epsilon = 0.02
g = 5.0
k = 8.0
uth = 0.1

dx = 0.2
dy = 0.2
dt = 0.01
L_x = 100
L_y = 100
frame_step = 100
                
reset_times = [2850, 2375, 1900, 1425, 950]
reset_heights =  [70, 60, 50, 40, 30]
               
# Modify electrode_positions to be an array of labels and a list of arrays for positions
electrode_labels = np.array(["L", "R", "V1", "V2", "V3", "V4", "V5", "V6", "LL"])
electrode_positions = np.array([
    [150, 130, 40],  # left shoulder
    [-150, 130, 40],  # right shoulder
    [-15, 15, 40],
    [15, 15, 40],
    [45, 0, 40],
    [80, -7, 40],
    [110, -15, 40],
    [130, -20, 40],
    [125, -150, 40]
])

# Initialize permanent variables outside the function
permanent_distance_arrays = []
permanent_ecg_arrays = []

#Set the number of cores to run on (only on own hardware)
#set_num_threads(3)

@njit
def calc_distances(L_x, L_y, dx, dy, sim_time, frame_step, electrode_positions):
    distance_arrays = []
    ecg_arrays = []

    for i in range(len(electrode_positions)):
        x, y, z = electrode_positions[i]
        distance = np.zeros((L_y, L_x))

        #array of zeros to store ecg values every recorded frame_step
        ecg_values = np.zeros((int((sim_time/frame_step))+1))

        for i in range(L_x):
            for j in range(L_y):
                distance[j, i] = np.sqrt((i - x)**2 + (j - y)**2 + (0 - z)**2)

        distance_arrays.append(distance)
        ecg_arrays.append(ecg_values)

    #return counter for number of frames passed
    #frame_count = 0

    return distance_arrays, ecg_arrays


@jit(parallel=True, nopython=True)
def simulate(time_run, reset_time, reset_height, electrode_positions):
    u = np.zeros((L_x, L_y))
    n = np.zeros((L_x, L_y))
    u_history = [u]
    u_next = np.zeros((L_x, L_y))

    u[:, :3] = 1
    electrode_distances, electrode_signals = calc_distances(L_x, L_y, dx, dy, time_run, frame_step, electrode_positions)
    frame_count = 0

    for t in range(time_run):
        if t == reset_time:
            print(f"Creating rotor at: {t}")
            u[(reset_height):, :] = 0
            u_history = [u]

        for i in prange(L_y -1):
            for j in range(L_x-1): 
                # Find central difference for each element
                CD_x = (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / dx ** 2
                CD_y = (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / dy ** 2

                # Update equations
                f = -k * u[i, j] * (u[i, j] - uth)  * (u[i, j] - 1)
                n[i, j] = n[i, j] + dt * (epsilon * (g * u[i,j] - n[i,j]))
                u_next[i, j] = u[i,j] + dt * (D * (CD_x + CD_y) + f - n[i,j])

                # boundary conditions
                u_next[0,:] = u_next[1,:]
                u_next[L_x-1,:] = u_next[L_x-2,:]
                u_next[:,0] = u_next[:,1]
                u_next[:,L_x-1] = u_next[:,L_x-2]

                if t % frame_step == 0:
                    for idx in range(len(electrode_labels)):
                        integrand = (D * (CD_x + CD_y)) / electrode_distances[idx][i, j]
                        electrode_signals[idx][frame_count] += integrand



        u = np.copy(u_next)

        if t % frame_step == 0:
            frame_count += 1
            u_history.append(np.copy(u))
            print(f"Iteration {t}/{time_run}")

    return u_history, electrode_signals


def sim_all(reset_times, reset_heights, electrode_positions):
    current_dir = os.getcwd()

    for i in range(len(reset_times)):
        for j in range(len(reset_heights)):
            for k in range(-2, 2):
                for l in range(-2, 2):
                    sim_time = 4700 * 7 + reset_times[i]

                    reset = reset_times[i] + (k * 38)
                    height = reset_heights[j] + l

                    # Pass electrode_distances to simulate function
                    u_history, electrode = simulate( sim_time, reset, height, electrode_positions)

                    df_electrode_data = pd.DataFrame({label: data for label, data in zip(electrode_labels, electrode)})

                    if current_dir:
                        csv_filename = os.path.join(current_dir,
                                                   f'electrode_data_reset_{reset_times[i]}_height_{reset_heights[j]}_x_{k}_y_{l}.csv')
                        df_electrode_data.to_csv(csv_filename, index=False)
                        print(f"Electrode data saved to: {csv_filename}")


# Call sim_all function
sim_all(reset_times, reset_heights, electrode_positions)

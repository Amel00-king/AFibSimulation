'''
Running Fitz-Hugh Nagumo simulations on a 2D grid and the associated ECG signals 
'''
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd 
from google.colab import files



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


''' only returns 100 values per electrode --> attempted fix done'''
def calc_distances(electrodes, L_x, L_y, dx, dy):
    
    distances = {label: np.zeros((L_y, L_x)) for label in electrode_positions.keys()}
    
    #calculate distance from each electrode to every element on the 2D 100x100 simulation grid
    for label, pos in electrodes.items():
        x, y, z = pos
        for i in range(L_x):
            for j in range(L_y):

                distances[label][i,j] = np.sqrt((i - x)**2 + (j - y)**2 + (0 - z)**2)

    
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

        for i in range(L_x - 1):
            for j in range(L_y -1):
            
                #find central difference for each element
                CD_x = (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / dx**2
                CD_y = (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / dy**2

                #update equations 
                f = -k * u[i,j] * (u[i,j] - uth) * (u[i,j] - 1)
                n[i,j] = n[i,j] + dt * (eps * (g * u[i,j] - n[i,j]))
                u_next[i,j] = u[i,j] + dt * (D * (CD_x + CD_y) + f - n[i,j])

                # boundary conditions
                u_next[0,:] = u_next[1,:]
                u_next[L_x-1,:] = u_next[L_x-2,:]
                u_next[:,0] = u_next[:,1]
                u_next[:,L_x-1] = u_next[:,L_x-2]

                
                #itterate through electrodes
                for label, distance in electrode_distances.items():
                    x, y, z = electrodes[label]
                    
                    #calculate integrand 
                    integrand = (D * (CD_x + CD_y)) / distance

                    # Trapezoidal rule integration
                    integrated_signal = np.trapz(integrand)

                    # Save the integrated signal for the current timestep
                    electrode_data[label].append(integrated_signal)
        
        
        if t%100 == 0:
            #save frames of u to u_history to later create gif of simulation
            u_history.append(np.copy(u))
    #
    u = np.copy(u_next)
        
    #convert electrode data to a DataFrame
    df = pd.DataFrame(electrode_data)
    df.to_csv('test.csv', encoding = 'utf-8-sig') 
    files.download('output.csv')

    # Plot electrode positions
    plot_electrodes(electrodes)

    return u_history, df



def simulate_all(electrode_pos, reset_times, reset_heights):

    for i in range(len(reset_times)):
        #It takes approx 4500-4600 time steps for a rotor wave to complete 1 cycle
        # therefore, total simulation time is 4700* 5 completed rotors, plus whatever time step the wave was initiated at
        sim_time = 4700 * 5 + reset_times(i)
        simulate(sim_time, reset_times(i), reset_heights(i), electrodes = electrode_pos)

    return  


# Create a function to update the plot at each frame
def update(frame):
    plt.clf()  # Clear the previous frame
    plt.imshow(u_history[frame], cmap='viridis')
    plt.title(f'Frame {frame}')


def visualize(u_history, electrode_df):
    # Create the figure and initial plot
    fig, ax = plt.subplots()
    im = ax.imshow(u_history[0], cmap='viridis')
    plt.title('Frame 0')

    # Create the animation
    num_frames = len(u_history)
    ani = FuncAnimation(fig, update, frames=num_frames, repeat=False)

    # Save the animation as a GIF
    ani.save('rotor.gif', writer='pillow', fps=8)  # Adjust fps as needed



    #plot all electrodes on same fig
    plt.figure()
    for label in df.columns:
        plt.plot(df.index, df[label], label=label)

    plt.xlabel("Time Step")
    plt.ylabel("Integrated Signal")
    plt.title("ECG Signals Over Time - All Electrodes")
    plt.legend()
    plt.show()


''' #check if the calc_distances is working for dictionaries
dist = calc_distances(electrode_positions,L_x,L_y,dx,dy)
# Assuming distances is the dictionary returned by calc_distances
for label, distance_array in dist.items():
    num_values = distance_array.size
    print(f"The number of values for label {label} is: {num_values}")
    num_rows, num_columns = distance_array.shape
    print(f"For label {label}, the array has {num_rows} rows and {num_columns} columns.")

'''



u_history, df = simulate(1050, reset_times[0], reset_heights[0],electrode_positions)
visualize(u_history,df) 
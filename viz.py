import scipy.io
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the matrix from the .mat file
data = scipy.io.loadmat('eptopic_r50_c50.mat')['matrix_to_export']

# Create a function to update the plot in the animation
def update(frame):
    plt.clf()
    plt.imshow(data[:, :, frame], cmap='viridis') # Adjust the colormap as needed

# Create an animation
ani = FuncAnimation(plt.gcf(), update, frames=data.shape[2], interval=100) # Adjust the interval as needed

# Save the animation as a GIF
ani.save('output.gif', writer='pillow') # Requires the Pillow library

plt.show()
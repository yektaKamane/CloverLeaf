import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def get_unique_j_values(file_pattern):
    j_values = set()
    directory = '../output/'
    for filename in os.listdir(directory):
        if filename.startswith(file_pattern) and filename.endswith('.csv'):
            parts = filename.split('_')
            j = int(parts[3].split('.')[0])  # Extract the j value
            j_values.add(j)
    return sorted(j_values)

def find_highest_j(file_pattern, i):
    j_values = set()
    directory = '../output/'
    for filename in os.listdir(directory):
        if filename.startswith(f'{file_pattern}_{i}_') and filename.endswith('.csv'):
            parts = filename.split('_')
            j = int(parts[3].split('.')[0])  # Extract the j value
            j_values.add(j)
    return max(j_values) if j_values else None

def create_grid(j, file_pattern):
    arrays = []
    directory = '../output/'
    i_values = sorted(set(int(filename.split('_')[2]) for filename in os.listdir(directory) if filename.startswith(file_pattern) and filename.endswith(f'_{j}.csv')))
    
    for i in i_values:
        file_name = f'{directory}/{file_pattern}_{i}_{j}.csv'
        if not os.path.exists(file_name):
            highest_j = find_highest_j(file_pattern, i)
            if highest_j is not None:
                file_name = f'{directory}/{file_pattern}_{i}_{highest_j}.csv'
        if os.path.exists(file_name):
            array = pd.read_csv(file_name).values
            arrays.append(array)
        else:
            # Handle the case where file doesn't exist or couldn't be found
            print(f"File {file_name} not found or does not exist.")
    
    # Ensure arrays has exactly 4 elements, using None for missing files
    while len(arrays) < 4:
        highest_j = find_highest_j(file_pattern, i_values[-1])  # Use highest j available for last i
        if highest_j is not None:
            file_name = f'{directory}/{file_pattern}_{i_values[-1]}_{highest_j}.csv'
            if os.path.exists(file_name):
                array = pd.read_csv(file_name).values
                arrays.append(array)
            else:
                arrays.append(None)
        else:
            arrays.append(None)
    
    # Filter out None values and concatenate existing arrays
    arrays = [arr for arr in arrays if arr is not None]
    
    if len(arrays) < 4:
        return None
    
    left_column = np.vstack((arrays[0], arrays[1]))
    right_column = np.vstack((arrays[2], arrays[3]))
    grid = np.hstack((left_column, right_column))
    
    return grid

def gather_grids(file_pattern):
    grids = []
    j_values = get_unique_j_values(file_pattern)
    
    for j in j_values:
        grid = create_grid(j, file_pattern)
        if grid is not None:
            grids.append(grid)
    
    return grids

def visualize_heatmap_anim(ax, array, vmin=0, vmax=3):
    cax = ax.imshow(array, cmap='coolwarm', aspect='auto', vmin=vmin, vmax=vmax)
    ax.axhline(y=array.shape[0] // 2, color='white', linestyle='--')
    ax.axvline(x=array.shape[1] // 2, color='white', linestyle='--')
    return cax

def animate_heatmaps(arrays, interval=500, save_as_gif=False, gif_name='heatmap_animation.gif'):
    fig, ax = plt.subplots()
    
    def update(frame):
        ax.clear()
        cax = visualize_heatmap_anim(ax, arrays[frame])
        ax.set_title('Heat Map of the 2D Array')
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
    
    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(arrays), interval=interval, repeat=False)
    
    if save_as_gif:
        writer = PillowWriter(fps=1000 / interval)
        ani.save(gif_name, writer=writer)
    else:
        plt.show()
    
    return ani

if __name__ == "__main__":
    file_pattern = 'array_data'
    grids = gather_grids(file_pattern)
    animate_heatmaps(grids, interval=50, save_as_gif=True, gif_name='heatmap_new.gif')

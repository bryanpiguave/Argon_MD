from RDF import Trajectory
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

def update(frame, ax, atom_trajectory):
    ax.clear()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Atom Trajectory with Trace - Frame: {frame}')

    # Plot the trace of the atom's path up to the current frame
    trace_x = atom_trajectory[:frame+1, 0]
    trace_y = atom_trajectory[:frame+1, 1]
    trace_z = atom_trajectory[:frame+1, 2]
    ax.plot(trace_x, trace_y, trace_z, color='gray', linewidth=1, linestyle='-')

    # Plot the current position of the atom with a marker
    current_position = atom_trajectory[frame]
    x_current = current_position[0]
    y_current = current_position[1]
    z_current = current_position[2]
    ax.scatter(x_current, y_current, z_current, s=50, c='blue', marker='o')

    # The axis limits are still set in the main function


def main():
    trajectory_file = '/home/bryan/Molecular_Dynamics/Project2/argon1_t_075.xyz'
    traj = Trajectory(trajectory_file, last_n=None, resolution=400) # Load all frames

    atom_index = 1  # Choose the index of the atom you want to track (0-based)
    atom_trajectory = traj.coordinates[:, atom_index, :]

    # Plotting each axis    
    plt.figure(figsize=(10, 6))
    plt.plot(atom_trajectory[:, 0], atom_trajectory[:, 1], label='X vs Y', color='blue')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('X vs Y Trajectory')
    plt.savefig('x_vs_y_trajectory.png')



    # # Determine axis limits based on the single atom's trajectory
    # x_min, x_max = np.min(atom_trajectory[:, 0]), np.max(atom_trajectory[:, 0])
    # y_min, y_max = np.min(atom_trajectory[:, 1]), np.max(atom_trajectory[:, 1])
    # z_min, z_max = np.min(atom_trajectory[:, 2]), np.max(atom_trajectory[:, 2])

    # # Add a small buffer to the limits for better visualization
    # buffer = 1.0
    # ax_xlim = [x_min - buffer, x_max + buffer]
    # ax_ylim = [y_min - buffer, y_max + buffer]
    # ax_zlim = [z_min - buffer, z_max + buffer]

    # # Animation of the trajectory in 3D
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.set_xlim(ax_xlim)
    # ax.set_ylim(ax_ylim)
    # ax.set_zlim(ax_zlim)

    # ani = FuncAnimation(fig, update, frames=traj.n_steps, fargs=(ax, atom_trajectory), interval=200)
    # # Save the animation
    # ani.save('single_atom_trajectory.gif', writer='imagemagick', fps=30)
    return



if __name__ == "__main__":
    main()
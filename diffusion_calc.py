import numpy as np
from RDF import Trajectory
import matplotlib.pyplot as plt
import argparse
import os
from scipy.stats import linregress
"""
Pick the coordinate of a certain atom and compute the mean squared displacement as a
#function of time. Then, calculate the diffusion constant for different temperatures and
#densities – you choose what temperature and density to use. Compare results
"""
parser = argparse.ArgumentParser(description='Calculate diffusion constant from MSD.')
parser.add_argument('--filename', type=str, help='Path to the trajectory file', default='/home/bryan/Molecular_Dynamics/Project2/argon09_t_075.xyz')
parser.add_argument('--skip', type=int, default=100, help='Number of snapshots to skip')
args = parser.parse_args()

def calculate_diffusion_coefficient(position_vector, dt):
    """
    Calculates the diffusion coefficient from a single particle trajectory.

    Args:
        position_vector (np.array): A numpy array of shape (n_steps, 3)
                                     representing the position of the particle over time.
        dt (float): The time step between each recorded position.

    Returns:
        float: The estimated diffusion coefficient.
    """
    n_steps = position_vector.shape[0]
    if n_steps < 2:
        return 0.0  # Cannot calculate MSD with less than 2 points

    initial_position = position_vector[0]
    squared_displacements = np.sum((position_vector - initial_position)**2, axis=1)
    time_array = np.arange(n_steps) * dt

    # We consider the MSD from the second point onwards (t > 0)
    time_array = time_array[1:]
    msd = squared_displacements[1:]

    if len(time_array) < 2:
        return 0.0

    # Perform a linear fit of MSD vs time
    slope, intercept, r_value, p_value, std_err = linregress(time_array, msd)
    diffusion_coefficient = slope / 6.0

    # Optional: Plot MSD vs Time
    plt.figure(figsize=(8, 6))
    plt.plot(time_array, msd, label='MSD')
    plt.plot(time_array, slope * time_array + intercept, 'r--', label=f'Linear fit (slope = {slope:.4f})')
    plt.xlabel('Time')
    plt.ylabel('Mean Squared Displacement')
    plt.title('Mean Squared Displacement vs Time')
    plt.legend()
    plt.grid(True)
    plt.show()

    return diffusion_coefficient


def main():
    # Load trajectory
    traj = Trajectory(args.filename, args.skip)
    print(traj.coordinates.shape)
    # Select the first atom
    atom_index = 0
    atom_positions = traj.coordinates[:, atom_index, :]
    # Calculate the time step
    dt = 0.00005  # Example time step in ps
    # Calculate the diffusion coefficient
    diffusion_constant = calculate_diffusion_coefficient(atom_positions, dt)

    # Print the diffusion constant in angstrom^2/ps
    print(f"Diffusion constant: {diffusion_constant:.4f} Å^2/ps")

if __name__ == "__main__":
    main()
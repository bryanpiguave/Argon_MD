import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm  # For progress bars
from scipy.spatial.distance import pdist
from collections import deque
import pandas as pd

from plot_aesthetics import axis_fontdict


FILENAME = 'argon1.xyz'
LAMMPS_RDF = '/home/bryan/Molecular_Dynamics/Project2/argon_rdf.dat'
parser = argparse.ArgumentParser(description='Compute and plot RDF from trajectory file.')
parser.add_argument('--filename', type=str, help='Path to the trajectory file', default=FILENAME)
parser.add_argument('--last_n', type=int, default=1, help='Number of last snapshots to keep')
parser.add_argument('--resolution', type=int, default=1000, help='Number of points in the final RDF')
parser.add_argument('--box_size', nargs=3, type=float, default=(5.39, 5.39, 5.39), help='Box size (x, y, z)')
parser.add_argument('--output_dir', type=str, default='.', help='Directory to save the plots')
parser.add_argument('--lammps_rdf', type=str, default=LAMMPS_RDF, help='Path to LAMMPS RDF file for comparison')
args = parser.parse_args()

class Trajectory:
    def __init__(self, filename, last_n=20, resolution=700):
        with open(filename, 'r') as f:
            data = f.readlines()
        self.n_atoms = int(data[0].split()[0])
        self.n_steps_total = int(len(data) / (self.n_atoms + 2))
        self.atom_list = [line.split()[0] for line in data[2:self.n_atoms+2]]
        self.last_n = last_n
        self.coordinates_queue = deque(maxlen=last_n)
        self.resolution = resolution
        
        # Read all steps but only keep the last_n in the queue
        for step in range(self.n_steps_total):
            i = step * (self.n_atoms + 2)
            snapshot = np.array([
                [float(val) for val in line.split()[1:4]] 
                for line in data[i+2:i+self.n_atoms+2]
            ])
            self.coordinates_queue.append(snapshot)
        
        # Convert the queue to a numpy array for easier processing
        self.coordinates = np.array(self.coordinates_queue)
        self.n_steps = len(self.coordinates)

    def compute_distance(self, a, b, box_size):
        """Minimum image convention distance calculation"""
        delta = np.abs(a - b)
        delta = np.where(delta > 0.5 * box_size, delta - box_size, delta)
        return np.sqrt(np.sum(delta**2))

    def compute_rdf(self, box_size=None):
        box_size = np.array(box_size)
        r_cutoff = 8 
        dr = r_cutoff / self.resolution
        hist = np.zeros(self.resolution)
        
        avg_density = self.n_atoms / np.prod(box_size)
        
        for i in tqdm(range(self.n_steps), desc="Computing RDF", unit="step"):
            snapshot = self.coordinates[i]
            # Calculate pairwise distances
            dists = pdist(snapshot, metric='euclidean')
            # Bin the distances
            hist += np.histogram(dists, bins=self.resolution, range=(0, r_cutoff))[0]
        
        # Normalize the histogram to get g(r)
        # Volume of each shell
        shell_volume = 4.0 / 3.0 * np.pi * ((np.arange(1, self.resolution + 1) * dr)**3 - ((np.arange(self.resolution) * dr)**3))
        # Correct for the volume of the shell
        shell_volume = np.where(shell_volume > 0, shell_volume, 1)
        # Normalize by the average density and the volume of the shell
        hist = hist / (self.n_steps*shell_volume )
        # Moving average
        hist = np.convolve(hist, np.ones(25)/25, mode='same')
        
        
        # Normalize by the number of atoms in the shell
        self.radii = np.arange(0, r_cutoff, dr)
        self.g_of_r = hist


    def plot_rdf(self, filename=""):
        plt.figure(figsize=(10, 6))
        plt.plot(self.radii, self.g_of_r, linewidth=2)
        plt.xlabel('r (Å)', fontsize=14)
        plt.ylabel('g(r)', fontsize=14)
        plt.title('Radial Distribution Function', fontsize=16)
        plt.grid(alpha=0.3)
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')

def main():
    trajectory = Trajectory(args.filename, args.last_n, args.resolution)
    print(f"Loaded {trajectory.n_steps} steps with {trajectory.n_atoms} atoms each.")
    print("Number of total steps:", trajectory.n_steps_total)
    print(trajectory.coordinates.shape)
    
    trajectory.compute_rdf(args.box_size)
    trajectory.plot_rdf(f"{args.output_dir}/rdf_plot.png")

    if args.lammps_rdf:
        with open(args.lammps_rdf, 'r') as f:
            lines = f.readlines()
        # Skip the first 4 lines of the LAMMPS RDF file
        lines = lines[4:]    
       # Read the LAMMPS RDF file in chunks
        with open(args.lammps_rdf, 'r') as f:
            lines = f.readlines()
            lines = lines[4:]  # Skip the first 4 lines
            df = pd.DataFrame(columns=['step','r', 'g_r', 'g_r_err'])
            df = pd.read_csv(args.lammps_rdf, delim_whitespace=True, header=None, names=['step','r', 'g_r', 'g_r_err'], skiprows=4)
            # Convert the columns to numeric, forcing errors to NaN
            df['step'] = pd.to_numeric(df['step'], errors='coerce')
            df['r'] = pd.to_numeric(df['r'], errors='coerce')
            df['g_r'] = pd.to_numeric(df['g_r'], errors='coerce')   
            df['g_r_err'] = pd.to_numeric(df['g_r_err'], errors='coerce')
        #Drop where step is 320200
        print("Dropping step 320200 from LAMMPS RDF data")
        df = df.drop(df[df['r'] == args.resolution].index)
        # Group by 'r' and calculate the mean and standard deviation

        avg_df = df.groupby('r').agg({'g_r': 'mean', 'g_r_err': 'std'}).reset_index()
        avg_df = avg_df.reset_index(drop=True)
        avg_df['r'] = avg_df['r'].astype(float)
        avg_df['g_r'] = avg_df['g_r'].astype(float)
        avg_df['g_r_err'] = avg_df['g_r_err'].astype(float)
        # Moving average
        window_size = 10
        avg_df['g_r'] = avg_df['g_r'].rolling(window=window_size, min_periods=1).mean()
        avg_df['g_r_err'] = avg_df['g_r_err'].rolling(window=window_size, min_periods=1).mean()

        # Plot the LAMMPS RDF
        plt.figure(figsize=(10, 6))
        # Plot the average g(r) 
        plt.plot(avg_df['r'], avg_df['g_r'], linewidth=2, label='LAMMPS RDF', color='darkblue')
        plt.plot(trajectory.radii, trajectory.g_of_r, linewidth=2, label='Computed RDF', color='orange')
        plt.xlabel('r (Å)', fontdict=axis_fontdict)
        plt.ylabel('g(r)', fontdict=axis_fontdict)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.title('Radial Distribution Function from LAMMPS', fontsize=24, fontweight='bold', color='darkblue')
        plt.grid(alpha=0.3)
        plt.legend(fontsize=16)
        plt.savefig(f"{args.output_dir}/lammps_rdf_plot.png", dpi=300, bbox_inches='tight')




    

if __name__ == "__main__":
    main()
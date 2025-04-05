import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm  # For progress bars
from scipy.spatial.distance import pdist
from collections import deque
import pandas as pd
FILENAME = 'argon1.xyz'
LAMMPS_RDF = '/home/bryan/Molecular_Dynamics/Project2/argon_rdf.dat'
parser = argparse.ArgumentParser(description='Compute and plot RDF from trajectory file.')
parser.add_argument('--filename', type=str, help='Path to the trajectory file', default=FILENAME)
parser.add_argument('--last_n', type=int, default=1, help='Number of last snapshots to keep')
parser.add_argument('--resolution', type=int, default=400, help='Number of points in the final RDF')
parser.add_argument('--box_size', nargs=3, type=float, default=(5.39, 5.39, 5.39), help='Box size (x, y, z)')
parser.add_argument('--output_dir', type=str, default='.', help='Directory to save the plots')
parser.add_argument('--lammps_rdf', type=str, default=LAMMPS_RDF, help='Path to LAMMPS RDF file for comparison')
args = parser.parse_args()

class Trajectory:
    def __init__(self, filename, last_n=10, resolution=400):
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
        r_cutoff = min(box_size) / 2.0
        
        dr = r_cutoff / self.resolution
        hist = np.zeros(self.resolution)
        
        # Precompute shell volumes and radii
        r = np.linspace(dr/2, r_cutoff-dr/2, self.resolution)
        shell_vol = 4 * np.pi * r**2 * dr
        
        if box_size is None:
            # For non-periodic systems, compute density based on the actual volume
            # (assuming coordinates are in a box of size max_coord - min_coord)
            min_coords = np.min(self.coordinates, axis=(0, 1))
            max_coords = np.max(self.coordinates, axis=(0, 1))
            volume = np.prod(max_coords - min_coords)
            avg_density = self.n_atoms / volume
        else:
            avg_density = self.n_atoms / np.prod(box_size)
        
        for step in range(self.n_steps):
            if box_size is None:
                dists = pdist(self.coordinates[step])
            else:
                dists = pdist(self.coordinates[step], 
                            lambda u, v: self.compute_distance(u, v, box_size))
            hist += np.histogram(dists, bins=self.resolution, range=(0, r_cutoff))[0]
        
        norm = shell_vol * avg_density * self.n_steps * self.n_atoms
        g_of_r = hist / norm
        
        self.radii = r
        self.g_of_r = g_of_r

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
    trajectory.compute_rdf(args.box_size)
    trajectory.plot_rdf(f"{args.output_dir}/rdf_plot.png")

    if args.lammps_rdf:
        # with open(args.lammps_rdf, 'r') as f:
        #     lines = f.readlines()

        # # Skip the first 4 lines of the LAMMPS RDF file
        # lines = lines[4:]
        # chunk_size = 500
        # list_of_dataframes = []
        
        # # Skip a line  for every 300 lines collected
        # counter = 0
        # for i in tqdm(range(0, len(lines)), desc="Processing LAMMPS RDF file"):
        #     if counter % chunk_size == 0:
        #         chunk = lines[i:i+chunk_size]
        #         # Convert the chunk to a DataFrame
        #         df = pd.DataFrame([line.split() for line in chunk])
        #         list_of_dataframes.append(df)
        #         # Skip the next line
        #         line = lines.readline()
                
        #     counter += 1
        # # for each dataframe get 
            

        # Load LAMMPS RDF data for comparison
        df = pd.read_csv(args.lammps_rdf,
                 skiprows=4,          # Skip 4 header lines
                 delim_whitespace=True, 
                 header=None)
        
        # Columns: Row, c_myRDF[1], c_myRDF[2], c_myRDF[3]
        r = df[1].values   # Distance (r)
        g_r = df[2].values # RDF (g(r))

        # Moving average
        window_size = 100
        moving_avg = pd.Series(g_r).rolling(window=window_size, min_periods=1).mean()
        g_r = moving_avg.values
        
        plt.figure(figsize=(10, 6))
        plt.plot(r, g_r, label='LAMMPS RDF', linewidth=2,alpha=0.8)
        
        plt.xlabel('r (Å)', fontsize=14)
        plt.ylabel('g(r)', fontsize=14)
        plt.title('Comparison of RDF', fontsize=16)
        plt.legend()
        plt.grid(alpha=0.3)
        plt.savefig(f"{args.output_dir}/rdf_comparison.png", dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()
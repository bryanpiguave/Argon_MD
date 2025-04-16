import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm  # For progress bars
import pandas as pd
import numpy as np
from movement import Trajectory 
from plot_aesthetics import axis_fontdict


LAMMPS_RDF = 'argon_rdf.dat'
parser = argparse.ArgumentParser(description='Compute and plot RDF from trajectory file.')
parser.add_argument('--filename', type=str, help='Path to the trajectory file')
parser.add_argument('--last_n', type=int, default=10000, help='Number of last snapshots to keep')
parser.add_argument('--resolution', type=int, default=6000, help='Number of points in the final RDF')
parser.add_argument('--box_size', nargs=3, type=float, help='Box size (x, y, z)')
parser.add_argument('--output_dir', type=str, default='.', help='Directory to save the plots')
parser.add_argument('--lammps_rdf', type=str, default=LAMMPS_RDF, help='Path to LAMMPS RDF file for comparison')
parser.add_argument('--reduced_density', type=float, default=0.5, help='Reduced density for RDF calculation')
args = parser.parse_args()




def main():
    trajectory = Trajectory(args.filename, args.resolution)
    print(f"Loaded {trajectory.n_steps} steps with {trajectory.n_atoms} atoms each.")
    print(trajectory.coordinates.shape)
    box_size = np.array(args.box_size)
    box_size = 4 * box_size
    trajectory.compute_rdf(box_size)
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
        radii = trajectory.radii
        g_of_r = trajectory.g_of_r
        # Moving average
        window_size = 10
        g_of_r = pd.Series(g_of_r).rolling(window=window_size, min_periods=1).mean()
        plt.plot(trajectory.radii, trajectory.g_of_r, linewidth=2, color='orange', alpha=0.3)
        plt.plot(radii  , g_of_r, linewidth=2, label='Computed RDF', color='orange', alpha=0.5)
        plt.plot(avg_df['r'], avg_df['g_r'], linewidth=2, label='LAMMPS RDF', color='darkblue')
        plt.xlabel('r (Å)', fontdict=axis_fontdict)
        plt.ylabel('g(r)', fontdict=axis_fontdict)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.title(r'Radial Distribution Function at ' +r'$\rho^* =$'+str(args.reduced_density) , fontsize=24, fontweight='bold', color='darkblue')
        plt.grid(alpha=0.3)
        plt.legend(fontsize=16)
        plt.savefig(f"{args.output_dir}/lammps_rdf_plot.png", dpi=300, bbox_inches='tight')


    

if __name__ == "__main__":
    main()

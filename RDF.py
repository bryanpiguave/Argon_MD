import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm  # For progress bars
import pandas as pd
import numpy as np
import sys # Used for progress indicator
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



class Trajectory:
    """
    Parses an entire XYZ trajectory file.

    Attributes:
        filename (str): Path to the XYZ file.
        n_atoms (int): Number of atoms per frame.
        atom_list (list): List of atom types/symbols from the first frame.
        coordinates (np.ndarray): A 3D NumPy array storing the coordinates
                                   of all steps.
                                   Shape: (n_steps, n_atoms, 3).
        n_steps (int): Total number of steps (frames) read and stored from the file.
        resolution (int): A parameter from the original class (usage not defined here).
                          Defaults to 700.
    """
    def __init__(self, filename, resolution=700):
        """
        Initializes the Trajectory object by reading the entire XYZ file.

        Args:
            filename (str): The path to the XYZ trajectory file.
            resolution (int, optional): A parameter for potential future use.
                                        Defaults to 700.
        """
        self.filename = filename
        self.resolution = resolution
        self.n_atoms = 0
        self.atom_list = []
        all_coordinates_list = []
        step_counter = 0 # Use this for counting steps read

        print(f"Reading trajectory file: {filename}...")

        try:
            with open(filename, 'r') as f:
                while True:
                    # --- Read Header ---
                    # 1. Read number of atoms (or check for EOF)
                    line_num_atoms = f.readline()
                    if not line_num_atoms:
                        break # End of file reached

                    try:
                        current_n_atoms = int(line_num_atoms.strip())
                    except ValueError:
                        print(f"Warning: Could not parse number of atoms at step {step_counter}. Line: '{line_num_atoms.strip()}'", file=sys.stderr)
                        # Attempt to skip this potentially corrupted frame
                        if self.n_atoms > 0: # If we know n_atoms from previous frames
                            for _ in range(self.n_atoms + 1): # Skip comment + atom lines
                                f.readline()
                        else: # If it's the first frame, we can't proceed
                           raise ValueError("Could not parse n_atoms on the first frame.")
                        continue # Skip to next potential frame header

                    if self.n_atoms == 0:
                        self.n_atoms = current_n_atoms # Set n_atoms based on the first frame
                    elif current_n_atoms != self.n_atoms:
                        print(f"Warning: Number of atoms changed at step {step_counter} ({current_n_atoms} vs {self.n_atoms}). Using first frame's count.", file=sys.stderr)
                        # Adapt logic if variable number of atoms is expected/handled differently

                    # 2. Read comment line (and ignore it for now)
                    f.readline() # Skip the comment/timestep line

                    # --- Read Coordinates for this step ---
                    # Ensure n_atoms has been set before allocating array
                    if self.n_atoms <= 0:
                         raise ValueError("Number of atoms is not positive, cannot read coordinates.")

                    coords_current_step = np.zeros((self.n_atoms, 3), dtype=float)
                    is_first_step = (step_counter == 0) # Check if it's the very first valid step

                    read_success = True
                    for i in range(self.n_atoms):
                        line_atom = f.readline()
                        if not line_atom:
                            print(f"Warning: Unexpected end of file while reading coordinates for step {step_counter}.", file=sys.stderr)
                            read_success = False
                            break # Exit the inner loop

                        parts = line_atom.split()
                        try:
                            if is_first_step:
                                self.atom_list.append(parts[0]) # Store atom type from first step
                            # Store coordinates
                            coords_current_step[i, 0] = float(parts[1]) # x
                            coords_current_step[i, 1] = float(parts[2]) # y
                            coords_current_step[i, 2] = float(parts[3]) # z
                        except (IndexError, ValueError) as e:
                            print(f"Warning: Error parsing line {i+3} of step {step_counter}: '{line_atom.strip()}'. Error: {e}", file=sys.stderr)
                            read_success = False
                            # Option: Fill with NaN or skip frame? For now, mark as failed.
                            coords_current_step[i, :] = np.nan
                             # If atom list wasn't filled yet during the first step, add a placeholder
                            if is_first_step and len(self.atom_list) <= i:
                                self.atom_list.append("?")


                    if read_success:
                        all_coordinates_list.append(coords_current_step)
                        step_counter += 1
                        # Simple progress indicator
                        if step_counter % 1000 == 0:
                            print(f"\rRead {step_counter} steps...", end="")

                print(f"\rRead {step_counter} steps... Done.")


        except FileNotFoundError:
            print(f"Error: File not found at {filename}", file=sys.stderr)
            self.coordinates = np.empty((0, 0, 3))
            self.n_steps = 0
            return # Stop initialization
        except Exception as e:
            print(f"An error occurred during file reading: {e}", file=sys.stderr)
            self.coordinates = np.empty((0, 0, 3))
            self.n_steps = 0
            return # Stop initialization


        # --- Convert the list of coordinates to the final NumPy array ---
        if not all_coordinates_list:
            print("Warning: No complete steps were read from the file.", file=sys.stderr)
            # Ensure n_atoms is somewhat valid before creating the empty array shape
            atom_dim = self.n_atoms if self.n_atoms > 0 else 0
            self.coordinates = np.empty((0, atom_dim, 3))
            self.n_steps = 0
        else:
             # Directly convert the whole list
             self.coordinates = np.array(all_coordinates_list)
             self.n_steps = len(self.coordinates) # n_steps is the total steps read
             print(f"Stored all {self.n_steps} steps found.")

    def compute_rdf(self, box_size):
        print('Box size:', box_size)
        box_size = np.array(box_size)
        r_cutoff = min(box_size) / 2.0
        dr = r_cutoff / self.resolution
        hist = np.zeros(self.resolution)
        
        avg_density = self.n_atoms / np.prod(box_size)
        
        for snapshot in tqdm(self.coordinates, desc="Computing RDF", unit="frame"):
            for i in range(self.n_atoms):
                for j in range(i + 1, self.n_atoms):
                    # Minimum image convention
                    delta = snapshot[j] - snapshot[i]
                    delta -= box_size * np.round(delta / box_size)
                    r = np.linalg.norm(delta)
                    
                    if r < r_cutoff:
                        bin_index = int(r / dr)
                        if bin_index < self.resolution:
                            hist[bin_index] += 2  # Each pair counts once for i-j and j-i
        
        radii = (np.arange(self.resolution) + 0.5) * dr
        shell_volumes = 4.0 * np.pi * radii**2 * dr
        N_ideal = avg_density * shell_volumes * self.n_atoms * self.n_steps
        g_r = hist / N_ideal
        
        self.radii = radii
        self.g_of_r = g_r


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

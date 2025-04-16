import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd 
from plot_aesthetics import axis_fontdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, default="output_files/grid_data/grid_data.txt")
parser.add_argument("--output_dir", type=str, default="output_files/grid_data")
args = parser.parse_args()
# Sample data with header
df = pd.read_csv(args.filename, sep='\s+', header=0)
#Drop duplicate rows
df = df.drop_duplicates()
def plot_diffusion_coefficients(df):
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot Diffusion vs Density at constant temperature
    temperature_filter = df['Temperature'] == 300
    df_filtered = df[temperature_filter]
    # Order by density
    df_filtered = df_filtered.sort_values(by='Density')
    ax1.plot(df_filtered['Density'], df_filtered['Diffusion_Coefficient'], 'ro-', linewidth=2, markersize=8)
    ax1.set_xlabel('Density (g/cm³)', fontdict=axis_fontdict)
    ax1.tick_params(axis='both', labelsize=18)
    ax1.set_ylabel(r'Diffusion Coefficient (×$10^{-4}$ cm²/s)', fontdict=axis_fontdict)
    ax1.set_title(r'$\mathbf{D}$  vs $\mathbf{\rho}$', fontsize=24, weight='bold', color='darkblue')
    ax1.grid(True, linestyle='--', alpha=0.6)
    
    # Plot Diffusion vs Temperature at density = 0.9
    density_filter = df['Density'] == 0.9
    df_filtered = df[density_filter]
    # Order by temperature
    df_filtered = df_filtered.sort_values(by='Temperature')
    ax2.plot(df_filtered['Temperature'], df_filtered['Diffusion_Coefficient'], 'go-', linewidth=2, markersize=8)
    ax2.set_xlabel('Temperature (K)',fontdict=axis_fontdict)
    ax2.tick_params(axis='both', labelsize=18)
    ax2.set_title(r'$\mathbf{D}$  vs $\mathbf{T}$', fontsize=24, weight='bold', color='darkblue')
    ax2.grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    plt.savefig(f"{args.output_dir}/diffusion_coefficients.png", dpi=300)

def main():
    plot_diffusion_coefficients(df)

if __name__ == "__main__":
    main()

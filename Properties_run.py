import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import os 
FILENAME = 'out.argon1'  # Default filename

parser = argparse.ArgumentParser(description='Process LAMMPS output data.')
parser.add_argument('--filename', type=str, default=FILENAME, help='Path to the LAMMPS output file')
parser.add_argument('--output_dir', type=str, default='.', help='Directory to save the plots')
args = parser.parse_args()
# Read file

with open(FILENAME, 'r') as f:
    data = f.readlines()

if os.path.exists(args.output_dir) == False:
    os.makedirs(args.output_dir,exist_ok=True)

data_chunks = []
current_chunk = []
in_data_section = False

for line in data:
    if 'Step          Temp' in line:
        if current_chunk:
            data_chunks.append(current_chunk)
        current_chunk = []
        in_data_section = True
    elif 'Loop' in line:
        if current_chunk:
            data_chunks.append(current_chunk)
        current_chunk = []
        in_data_section = False
    elif in_data_section and line.strip():
        current_chunk.append(line.split())

if current_chunk:
    data_chunks.append(current_chunk)

# Create DataFrames for each chunk
dataframes = {}

if data_chunks:
    column_names = ['Step', 'Temp', 'E_pair', 'E_mol', 'TotEng', 'Press'] #hardcoded column names.
    for i, chunk in enumerate(data_chunks):
        if chunk:
            try:
                df = pd.DataFrame(chunk, columns=column_names) #use hardcoded names.
                dataframes[f"stage_{i+1}"] = df
            except ValueError:
                print(f"Warning: Chunk {i} has incorrect number of columns. Skipping.")
        else:
            print(f"Warning: Chunk {i} is empty. Skipping.")

    # Example: Print the first few DataFrames for verification
    for key, df in list(dataframes.items())[:3]:
        print(f"{key}:\n{df}\n")
else:
    print("No data chunks found.")


def Plotting_variable(dataframes, variable_name, keyword, output_dir): 
    axis_fontdict = {
        'family': 'sans-serif',  # Or any other font family
        'color':  'darkblue',
        'weight': 'bold',
        'size': 18,
    }
    # Plotting the data
    # Plotf {keyword} vs Step for each DataFrame with moving average
    plt.figure(figsize=(10, 6))
    window_size = 100  # Adjust window size as needed
    for key, df in dataframes.items():
        try:
            df['Step'] = df['Step'].astype(int)
            df[variable_name] = df[variable_name].astype(float)

            # Calculate moving average
            moving_avg = df[variable_name].rolling(window=window_size, min_periods=1).mean()
            plt.plot(df['Step'], df[variable_name], color='blue', alpha=0.5)
            plt.plot(df['Step'], moving_avg,  color='blue')

        except KeyError:
            print(f"Warning: 'Step' or variable_name not found in {key}. Skipping.")
    
    plt.xlabel('Step',fontdict=axis_fontdict)
    plt.ylabel(f'{keyword}',fontdict=axis_fontdict)
    plt.title(f'{keyword} vs Step',fontsize=24, fontweight='bold',color='darkblue')
    plt.axvline(x=10000,label='Equilibration NVE',color='black',linestyle='--',linewidth=2)
    plt.axvline(x=20000,label='Production Run',color='green',linestyle='--',linewidth=2)
    plt.legend(fontsize=16,loc='lower right')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{keyword}_vs_step.png')


# Temperature
variable_name = 'Temp'
keyword = 'Temperature'
Plotting_variable(dataframes, variable_name, keyword, args.output_dir)
# Pressure
variable_name = 'Press'
keyword = 'Pressure'
Plotting_variable(dataframes, variable_name, keyword, args.output_dir)
# Potential Energy
variable_name = 'E_pair'
keyword = 'Potential Energy'
Plotting_variable(dataframes, variable_name, keyword, args.output_dir)
# Total Energy
variable_name = 'TotEng'
keyword = 'Total Energy'
Plotting_variable(dataframes, variable_name, keyword, args.output_dir)

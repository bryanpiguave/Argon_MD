# Argon Simulation using LAMMPS <img src="argon_gen.jpg" alt="Fancy Argon Atom" width="50" align="center">

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
  - [Calculating Radial Distribution Function (g(r))](#calculating-radial-distribution-function-gr)
  - [Comparing g(r) with LAMMPS Built-in](#comparing-gr-with-lammps-built-in)
  - [Calculating Mean Squared Displacement (MSD) and Diffusion Constant](#calculating-mean-squared-displacement-msd-and-diffusion-constant)
- [License](#license)
- [Author](#author)

## Introduction

This project involves performing molecular dynamics simulations of Argon using the Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS). The primary goals are to calculate and analyze the structural and dynamical properties of Argon under different thermodynamic conditions. Specifically, this project focuses on:

1.  Calculating the radial distribution function (g(r)) from atom coordinates obtained from LAMMPS simulations and comparing it with the built-in LAMMPS calculation.
2.  Computing the mean squared displacement (MSD) of Argon atoms and subsequently determining the diffusion constant at various temperatures and densities.

## Features

This project provides the following functionalities:

* LAMMPS input scripts for simulating Argon at different densities and temperatures using the Lennard-Jones potential.
* Python scripts to:
    * Parse LAMMPS dump files and calculate the radial distribution function (g(r)).
    * Plot the calculated g(r)  at a fixed temperature and fixed density.
    * Compare the calculated g(r) with the output from LAMMPS' built-in RDF computation.
    * Analyze LAMMPS output to calculate the mean squared displacement (MSD) of selected atoms as a function of time.
    * Calculate the diffusion constant from the linear regime of the MSD plot.
* Flexibility to choose different temperatures and densities for analysis.

## Getting Started

### Prerequisites

Before running the simulations and analysis, ensure you have the following software installed:

* **LAMMPS:** You need a working installation of LAMMPS on your system. You can find installation instructions on the official LAMMPS website: [https://lammps.sandia.gov/](https://lammps.sandia.gov/)
* **Python 3:** Python 3 is required for the analysis scripts.
* **NumPy:** A fundamental package for numerical computation in Python. Install it using `pip install numpy`.
* **Matplotlib:** A comprehensive library for creating static, interactive, and animated visualizations in Python. Install it using `pip install matplotlib`.

### Installation

1.  **Clone the repository (if applicable):** If you have the project files in a repository, clone it to your local machine.
2.  **Navigate to the project directory:** Open your terminal or command prompt and navigate to the directory containing the LAMMPS input scripts and Python analysis codes.

## Usage

The project is structured into three main tasks, each with corresponding LAMMPS input scripts and Python analysis scripts.

### Calculating Radial Distribution Function (g(r))

1.  **Run the LAMMPS simulation:** Execute the `in.argon1` script using LAMMPS from your terminal:
    ```bash
    lmp < in.argon1
    ```
    This script will perform short simulations at densities 0.9, 1.0, and 1.1 (at T=0.75) and output atom coordinates to `dump.coord.0.9`, `dump.coord.1.0`, and `dump.coord.1.1`.
2.  **Run the Python analysis script:** Execute the `calculate_gr.py` script:
    ```bash
    python RDF.py
    ```
    This script will read the coordinate dump files, calculate the radial distribution function for each density, and generate a plot named `gr_T0.75.png`.

### Comparing g(r) with LAMMPS Built-in

1.  **Run the LAMMPS simulation for RDF calculation:** Execute the `in.argon.rdf` script using LAMMPS:
    ```bash
    lmp < in.argon.rdf
    ```
    This script will perform simulations and use the `compute rdf` and `fix ave/spatial` commands to calculate and output the radial distribution function directly to files named `rdf.0.9.dat`, `rdf.1.0.dat`, and `rdf.1.1.dat`.
2.  **Run the Python comparison script:** Execute the `compare_gr.py` script:
    ```bash
    python compare_gr.py
    ```
    This script will plot the g(r) calculated by your `calculate_rdf` function and the g(r) read from the LAMMPS output files on the same axes for each density, generating a comparison plot named `compare_gr_T0.75.png`.

### Calculating Mean Squared Displacement (MSD) and Diffusion Constant

1.  **Modify the LAMMPS MSD script:** Open the `in.argon.msd` script and adjust the `variable rho equal ...` and `variable T equal ...` lines to set the desired density and temperature for your MSD calculation.
2.  **Run the LAMMPS simulation for MSD:** Execute the `in.argon.msd` script using LAMMPS:
    ```bash
    lmp < in.argon.msd
    ```
    This script will simulate the Argon system at your chosen conditions and output the mean squared displacement of the first atom to the file `msd.dat`.
3.  **Analyze the MSD data (you will need to write a Python script for this):** You will need to create a new Python script (e.g., `analyze_msd.py`) to:
    * Read the `msd.dat` file.
    * Extract the time steps and MSD values.
    * Plot the MSD as a function of time.
    * Identify the linear regime at longer times.
    * Perform a linear fit to the linear regime of the MSD plot.
    * Calculate the diffusion constant (D) using the relationship:
        ```
        MSD(t) = 6Dt + constant
        ```
        where the slope of the linear fit is equal to `6D`.
    * You can then repeat this process for different temperatures and densities by modifying the `in.argon.msd` script and re-running the simulation.

## License

This project is licensed under the [Apache License 2.0](https://opensource.org/licenses/Apache-2.0).

# Author 
Bryan Piguave 
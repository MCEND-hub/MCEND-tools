# MCEND-tools

Tooling library for MCEND. Requirements depend on the individual scripts.

# `hhg_spec_gen`

The `hhg_spec_gen` is a small fortran library that can be used to generate high harmonic spectra. You may use a sliding window for the Fourier transform to identify time spans during which specific frequencies dominate (originate in). Compile using the provided make file.
The library takes an `expec.t` file as input and writes to `ft_omega.dat` as output.

# `proto_hdf5_writer`

The `proto_hdf5_writer` is a small Fortran helper routine that reads and writes hdf5 files.

# Python and bash scripts

The following Python scripts are provided for your convenience:
- `check_energy_conserved.py`: Checks if the total energy is conserved throughout the propagation (within a threshold, currently set as 1.0e-3, takes in `expec.t` as input file)
- `compare_expec_v2.py`: Compares if two `expec.t` files are the same within a given threshold (currently set as `1.0e-5`)
- `dipole-calc.py`: Calculates the dipole moment given the nuclear masses and internuclear distance from the `expec.t` file, as well as the electronic expectation value
- `plot-run-results.py`: Plots the MCEND output from the `expec.t` file

The following Jupyter notebooks are provided for your convenience:
- `spectra-mcend.ipynb`: Construct and plot high-harmonic and absorption spectra from the `expec.t` file

The following bash scripts are provided for your convenience:
- `clean.sh`: A script that removes all MCEND output in the current folder and the MCEND executable from `../bin`
- `collect-data.sh`: This script collects the MCEND output for you and moves it to a specificed folder. Can also trigger the calculation of spectra.
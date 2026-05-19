
# The code and plotting scripts to reproduce data and plots from the paper "Complex gas flows in magnetised protoplanetary disks promote the formation of dust traps at low fragmentation velocities"

# To download the data from zenodo:



# To run the simulations

## fidss simulation

make 
./2DMC setup_fidss.par

## fidmhd simulation

make
./2DMC setup_fidmhd.par

## mhdvfrag100 simulation

make
./2DMC setup_mhdvfrag100.par

# The plotting script to reproduce the figures from the paper can be found in scripts/paperplots.ipynb


#!/bin/bash
# This script is used by generate_plot_data.py to generate the plots for the figures in the paper: 
# 'An Efficient MCMC Approach to Energy Function Optimization in Protein Structure Prediction'

# To run an individual protein simulation, please use fold.sh

# $1 is seq.npz filepath
# $2 is fasta filepath
# $3 is storage file path for output
# $4 is the target energy, stops l-bfgs at target energy
# $5 is the filepath to the reference pdb structure to calculate rmsd for protein
# $6 is number of restarts
# $7 2 is gd, 3 is altmh

python ../foldingmc/RosettaTRMC.py $1 $2 $3 $4 $5 -r $6 -m $7 #--save_chk
echo "output at $3";


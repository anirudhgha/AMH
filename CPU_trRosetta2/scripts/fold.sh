#!/bin/bash

# Script can be used as follows 
# ./fold.sh [mode] [protein name] [target energy to stop simulation at] [reference pdb for RMSD]

# setting target energy to -100000000 will have it run for num_samples. Modes are described in the readme. 


if [ -z "$1" ]; then
    echo "[error] Need minimization method: 0,1,2 (gradient descent l-bgfs) or 3 (alternative metropolis-hastings) \nex. sh fold.sh # T0955";
    exit;
fi
echo "Set method to : $1";

if [ -z "$2" ]; then
    echo "[error] Need to pass the protein name corresponding to folder name in ../proteins. Ensure folder contains seq.npz and fasta.fasta. \nex. $ ./fold.sh 3 T0955";
    exit;
fi

mode=$1 
protein=$2          # T0955
t_energy=$3         # -100000 (primarily used to run simulations until the energy matches that of another optimization method)
ref_pdb="../proteins/$protein/${protein}_ref.pdb"          # ../proteins/$protein/protein-name_ref.pdb
out_path="../out"   # 
 
echo "processing protein: $protein";
echo "mode: $mode";
rm ../proteins/$protein/$mode_out.
python ../folding_mc/RosettaTRMC.py ../proteins/$protein/seq.npz ../proteins/$protein/fasta.fasta $out_path -r 1 -m $mode --no-fastrelax -ref $ref_pdb -target_energy -1000000000
echo "output at ../proteins/$protein/$mode_out.pdb";
echo -en "\007"
exit;
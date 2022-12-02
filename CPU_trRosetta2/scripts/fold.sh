#!/bin/bash

if [ -z "$1" ]; then
    echo "[error] Need minimization method: gd (gradient descent l-bgfs) or mc (alternative metropolis-hastings) \nex. sh fold.sh amh T0955";
    exit;
fi
echo "Set method to : $1";

if [ -z "$2" ]; then
    echo "[error] Need to pass the protein name corresponding to folder name in ../proteins. Ensure folder contains seq.npz and fasta.fasta. \nex. $ sh foldmc.sh T0955";
    exit;
fi

mode=0;
if [ $1 != "gd" ]; then
    echo "it worked";
    mode=3;    
else
    mode=2;
fi

echo "processing protein: $2";
echo "mode: $mode";
rm ../proteins/$2/$1_out.pdb
python ../foldingmc/RosettaTRMC.py ../proteins/$2/seq.npz ../proteins/$2/fasta.fasta ../out -r 1 -m $mode
echo "output at ../proteins/$2/$1_out.pdb";
echo -en "\007"
exit;
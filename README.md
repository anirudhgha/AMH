# An Efficient MCMC Approach to Energy Function Optimization in Protein Structure Prediction

## Getting Started
This repo contains two codes: ```CPU_trRosetta2``` and ```GPU_AMH```. ```CPU_trRosetta2``` contains code that modifies [trRosetta2's code](https://github.com/RosettaCommons/trRosetta2) to run alternative optimization methods on CPU. Run ```scripts/fold.sh``` for an example script to run one of the below optimization schemes on a protein. 

For example, to run l-bfgs on T0955, run ```./fold.sh 2 T0955``` in the scripts folder. Replace '2' according to the mapping below to specify other optimization schemes.  

- mode = 0,1,2 - default lbfgs (2 is used in paper)
- mode = 3 - alternating metropolis hastings
- mode = 4 - lbfgs with localized mh restarts
- mode = 5 - traditional with annealing
- mode = 6 - MALA with L-bfgs w/ noisy restarts
- mode = 7 - traditional with noisy restarts

*Notes* This code was originally run in Unix (Linux-Ubuntu 18.04). If encountering a bad interpreter error, try running ```sed -i -e 's/\r$//' fold.sh```. 

Each method requires a set of proteins to optimize, their target energies (may be set to -10000000000 if you want to run for some number of samples), an output folder path, and a pdb reference to calculate the rmsd. RMSD is scored as CA-CA distances, rather than full atom rmsd. All data for a run is stored in the info.csv file output at the output folder location. 

The second major code, ```GPU_AMH``` is for GPU acceleration of the Alternating Metropolis Hastings algorithm. This was run on google colab: The same folder may be loaded to google colab to execute this code. The user will have to modify the paths in the *user parameter* section of gpuinf.ipynb. The code used to parse the GPU output is also included, though again the path variables may change based on the user's google drive file structure. 

The GPU acceleration requires providing as input a CSV file with information pertaining to the protein to be folded. The proteins used in the manuscript are provided, and additional proteins may be formatted by the user. 

The CSV must contain the following information: 
- predicted distogram
- initial state (x,y,z)
- indices of CA atoms in the protein (see order via pdb)
- indices of CB atoms in the protein 
- indices of C atoms in the protein 
- indices of N atoms in the protein 



## Inputs
To run simulations on a new protein, create a folder for the protein in the ```proteins``` folder. Place the corresponding `fasta.fasta` and `seq.npz` files (as those names) in that folder. To measure RMSD, you would need some experimentally derived pdb file for that protein. Then use fold.sh to link these elements together and get a folded structure in the ```out``` folder. 


## Outputs
### info.csv
Has all the information regarding the run, including the score history. All columns are repeated except for score history. They are all the final values at the end of the run. 

### energy plot
Energy plot is accurate for the localized mh regions. The flat regions are reflective of the samples generated during the gradient descent portion. However, we only have the final energy after the gradient descent minmover is complete. Therefore, the flat region is just the final energy repeated over the samples the gradient descent runs for. The GD does not actually instantly drop. 



## Bibtex citation
@article{ https://doi.org/10.48550/arxiv.2211.03193,\
  doi = {10.48550/ARXIV.2211.03193},\
  url = {https://arxiv.org/abs/2211.03193},\
  author = {Ghantasala, Lakshmi A. and Jaiswal, Risi and Datta, Supriyo},\
  keywords = {Biomolecules (q-bio.BM), Quantitative Methods (q-bio.QM), Computation (stat.CO), FOS: Biological sciences, FOS: Biological sciences, FOS: Computer and information sciences, FOS: Computer and information sciences},\
  title = {An Efficient MCMC Approach to Energy Function Optimization in Protein Structure Prediction},\
  publisher = {arXiv},\
  year = {2022},\
  copyright = {Creative Commons Attribution 4.0 International}}

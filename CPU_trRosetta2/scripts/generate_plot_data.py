
# HOW TO USE: run either gen_plot1() or gen_plot2() in the __main__ at the bottom of the code. 
# DESCRIPTION: Script for generating the data needed to make the plots in the paper. This will run the alternating-mh scheme on various 
# proteins and save the results in different folders. 


# python ../foldingmc/RosettaTRMC.py $1 $2 $3 -r $4 -m 3
# $1 is seq.npz filepath
# $2 is fasta filepath
# $3 is storage file path for output
# $4 is the number of restarts


import subprocess
import os


# plot multiple proteins: min energy vs min energy scatterplot, rmsd vs rmsd scatterplot, num samples per protein to min energy
def gen_plot():
    nruns_limit = 10                          # total runs 
    nrestarts = 100                     # noisy restarts within one run
    root_out_path = '../out'            # store output figures, plots, csv data
    proteins_to_skip = []               #['T0955', 'T0980s2', 'T0953s1', 'T1006','0']
    proteins = [    
                    # '../proteins/T0955', 
                    # '../proteins/T1006',
                    # '../proteins/T0884', 
                    # '../proteins/T0953s1', 
                    # '../proteins/T0974s1'
                    '../proteins/T0980s2', 
                    # '../proteins/T1008', 
                    # '../proteins/T1046s1', 
                    # '../proteins/T1059', 
                    # '../proteins/T1072s2'
                ]
    ref_structures = [    
                    # '../proteins/T0955/T0955.pdb', 
                    # '../proteins/T1006/6qek.pdb',
                    # '../proteins/T0884/5t87.pdb', 
                    # '../proteins/T0953s1/6f45.pdb', 
                    # '../proteins/T0974s1/6tri.pdb'
                    '../proteins/T0955/T0955.pdb', 
                    # '../proteins/T1008/6msp.pdb', 
                    # '../proteins/T1046s1/6px4.pdb', 
                    # '../proteins/T1059/T0955.pdb', 
                    # '../proteins/T1072s2/6hk8.pdb'
                ]
    nruns_start_per_protein = [0]*len(proteins)

    target_energies = [100000000]*len(proteins)
    
    # 2 is gd w/ nr
    # 3 is MH w/ locmh (altmh)
    # 4 is gd w/ locmh
    # 5 is MH w/annealing
    # 6 is gd w/ mhcriteron (MALA/HMC) w/ nr
    # 7 is MH w/annealing w/ nr 
    args_mode = 5
    for i,pdir in enumerate(proteins):
        if pdir.split('/')[-1] in proteins_to_skip:
            continue
        for run in range(nruns_start_per_protein[i], nruns_limit):
            print('PROCESSING PROTEIN: {} run: {}'.format(pdir, run))
            print('PROTEIN: ', pdir.split('/')[-1])

            pre_num_path = root_out_path+'/'+pdir.split('/')[-1]
            if not os.path.isdir(pre_num_path):
                os.mkdir(pre_num_path)

            out_path = root_out_path+'/'+pdir.split('/')[-1]+'/'+str(run)
            print(out_path)

            if not os.path.isdir(out_path):
                os.mkdir(out_path)
            fasta_filepath = pdir+'/fasta.fasta'
            npz_filepath = pdir+'/seq.npz'
            
            print('OUTPUT PATH FOR RESULTS: ', out_path)
            subprocess.check_call(f"sh fold_subscript.sh {npz_filepath} {fasta_filepath} {out_path} {target_energies[i]} {ref_structures[i]} {nrestarts} {args_mode}", shell=True)
            print('COMPLETED RUN {}'.format(run))



if __name__ == "__main__":
    gen_plot()

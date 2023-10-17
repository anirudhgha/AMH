# Arranges all the data into a single dict, and writes it to a csv file to be used as input to the GPU. 

import time
import numpy as np
import mdtraj as md
import pickle
import subprocess
import jdata as jd
import csv
from scipy.io import savemat
import matplotlib.pyplot as plt
from fxpmath import Fxp

from pseudo_fpga import fpga

# ###########################################################################################
# general MCMC Parameters
noise_var_range = [10, 1]
cycles = 40
segment_sizes = [10, 40, 82]
iterations = [20, 25, 150]
beta = 1
n_runs = 10 # number of repetitions of experiment

# protein Specific Details
dist_filename = r'../T0955/contacts/T0955.pickle' 
initial_state_pdb_filename = r'../T0955/T0955_unfolded.pdb'

# misc
data_filename = 'T0955_cyc' + str(cycles) + '_ite' + str(iterations) + '_seg' +str(segment_sizes) + '_data.json'
mat_filename = 'T0955_cyc' + str(cycles) + '_ite' + str(iterations) + '_seg' +str(segment_sizes) + '_data.mat'

server_ipv4_address = 'ubuntu@'
##############################################################################################


# All values to be sent to GPU stored in r
r={}

with open(dist_filename, 'rb') as infile:
    new_dict = pickle.load(infile, encoding='latin1')
r['distogram'] = new_dict['probs']
# else:
    # r['distogram'] = new_dict['probs'][dist_start - 1:dist_end, dist_start - 1:dist_end]

default_traj = md.load(initial_state_pdb_filename)
topology = default_traj.topology
r['topology'] = topology
table, bonds = topology.to_dataframe()
index = table.index
r['molecule_xyz'] = np.squeeze(np.array(default_traj.xyz[0, :, :]))
r['n_residues'] = default_traj.n_residues
r['n_atoms'] = default_traj.n_atoms
r['CA'] = np.array(index[table['name'] == 'CA'])  # folding
r['CB'] =  np.array(index[table['name'] == 'CB'])  # energy calculation
r['N'] =  np.array(index[table['name'] == 'N'])
r['C'] =  np.array(index[table['name'] == 'C'])
r['O'] =  np.array(index[table['name'] == 'O'])  # phi
r['H'] =  np.array(index[table['name'] == 'H'])  # psi
r['noise_var_range'] = np.array([10,1])
r['beta_cycles'] = np.array([1, cycles])
r['beta'] = 1
r['cycles']=cycles
r['segment_sizes'] = np.array(segment_sizes)
r['iterations'] = np.array(iterations)

jd.save(r, 'fpga_run_data/'+data_filename)
savemat(file_name='fpga_run_data/'+ mat_filename, mdict=r)


# flatten out the distogram into a 1d string to make reading easier in c++ (and sending to fpga easier). 
# I will have n_residue info to read it properly in the fpga

# uncomment below lines to get new data for GPU to run!


flattened_distogram = [] #np.zeros(1,r['n_residues']*r['n_residues'] * 64)
for i in range(r['n_residues']):
    for j in range(r['n_residues']):
        flattened_distogram += list(r['distogram'][i][j][:])

flattened_distogram = np.array(flattened_distogram)
print(f'flattened_distogram shape: {flattened_distogram.shape}')

flattened_log_distogram = np.log(flattened_distogram)
print(f"r[CA].shape: {r['CA'].shape}")

flattened_state = []
for mol in r['molecule_xyz']:
    for axis in mol:
        flattened_state += [axis]

print("initial state length: ", np.array(flattened_state).shape)



max_dist_weight = np.amax(flattened_distogram)
min_dist_weight = np.amin(flattened_distogram)

max_dist_logweight = np.amax(flattened_log_distogram)
min_dist_logweight = np.amin(flattened_log_distogram)

max_state = np.amax(flattened_state)
min_state = np.amin(flattened_state)

# print(max_dist_weight)
# print(min_dist_weight)
# print(max_dist_logweight)
# print(min_dist_logweight)
# print(max_state)
# print(min_state)


# convert ind array elements to 16 bit sequences
r['CA_bin'], r['CB_bin'], r['C_bin'], r['N_bin'], r['O_bin'] = [], [], [], [], []

for ind in r['CA']:
    binval = bin(ind)[2:]    
    binval = '0'*(16-len(binval)) + binval
    r['CA_bin'].append(binval)
for ind in r['CB']:
    binval = bin(ind)[2:]    
    binval = '0'*(16-len(binval)) + binval
    r['CB_bin'].append(binval)
for ind in r['C']:
    binval = bin(ind)[2:]    
    binval = '0'*(16-len(binval)) + binval
    r['C_bin'].append(binval)
for ind in r['N']:
    binval = bin(ind)[2:]    
    binval = '0'*(16-len(binval)) + binval
    r['N_bin'].append(binval)
for ind in r['O']:
    binval = bin(ind)[2:]    
    binval = '0'*(16-len(binval)) + binval
    r['O_bin'].append(binval)

# convert the dist and state to 32 bit fixed point sequences where [+-, 5 integer bits, 25 fractional bits]
flattened_distogram_bin, flattened_log_distogram_bin, flattened_state_bin  = [], [], []
for w in flattened_distogram:
    x = Fxp(w, signed=True, n_word=32, n_frac=25)
    flattened_distogram_bin.append(x.bin())
for w in flattened_log_distogram:
    x = Fxp(w, signed=True, n_word=32, n_frac=25)
    flattened_log_distogram_bin.append(x.bin())
for w in flattened_state:
    x = Fxp(w, signed=True, n_word=32, n_frac=25)
    flattened_state_bin.append(x.bin())

# write to csv file, one row per variable
f = open(data_filename, 'w')
writer = csv.writer(f)
writer.writerow(flattened_distogram_bin)
writer.writerow(flattened_log_distogram_bin)
writer.writerow(flattened_state_bin)
writer.writerow(r['CA_bin'])
writer.writerow(r['CB_bin'])
writer.writerow(r['C_bin'])
writer.writerow(r['N_bin'])
writer.writerow(r['O_bin'])
f.close()
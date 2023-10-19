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
noise_var_range = [10, 4]
cycles = 10
segment_sizes = [10000, 16]
iterations = [50, 50]
beta = 1
n_runs_per_protein = 1 # number of repetitions of experiment
n_proposal_limit = 100000 # if 0, run for cycles, else run until n_proposal_limit is reached. 

# 2653, 2847, 2689, 2859, 2763, 2859, 2712
# protein Specific Details (will search in file system, which is ../protein_name/etc...)
# T0955, T0953s2-D1, T0957s1-D2, T1019s1-D1, T0990-D1, T1006-D1, T0963-D2, T0991
protein_names = ['T0955', 'T0953s2', 'T0957s1', 'T1019s1', 'T0990', 'T1006','T0963', 'T0991', 'T0986s2']
dist_start = [1, 2, 38, 1,  1, 1, 40, 1, 1]
dist_end = [41, 45, 91, 58, 76, 77, 121, 118, 155]

protein_names = ['T0955']
dist_start = [1]
dist_end = [41]

strip_r_group = True
# misc
# data_filename = 'T0955_cyc' + str(cycles) + '_ite' + str(iterations) + '_seg' +str(segment_sizes) + '_fpgadata.json'
# mat_filename = 'T0955_cyc' + str(cycles) + '_ite' + str(iterations) + '_seg' +str(segment_sizes) + '_fpgadata.mat'

# server_ipv4_address = 'ubuntu@'
##############################################################################################

def strip_rgroup(xyz_coords):
    pass

for ind, name in enumerate(protein_names):
    dist_filename = r'../{}/contacts/{}.pickle'.format(name, name) #C:/Users/alasg/Documents/GitHub/Molecular_simulations/md/AlphaFold
    initial_state_pdb_filename = r'../{}/{}_unfolded.pdb'.format(name, name)
    
    # All values to be sent to fpga stored in r
    r={}

    with open(dist_filename, 'rb') as infile:
        new_dict = pickle.load(infile, encoding='latin1')
    # r['distogram'] = new_dict['probs']
    r['distogram'] = new_dict['probs'][dist_start[ind] - 1:dist_end[ind], dist_start[ind] - 1:dist_end[ind]]
    default_traj = md.load(initial_state_pdb_filename)
    
    # strip r group
    if strip_r_group:
        atom_names_to_keep = ['CA', 'CB', 'N', 'C']
        trunc_topology = default_traj.topology.copy()
        table, bonds = default_traj.topology.to_dataframe()
        deleted_atoms = []
        for i in range(default_traj.n_atoms - 1, -1, -1):
            if str(table['name'][i]) not in atom_names_to_keep:
                trunc_topology.delete_atom_by_index(i)
                deleted_atoms.append(i)
        table, bonds = trunc_topology.to_dataframe()
        index = table.index
        positions = np.squeeze(np.array(default_traj.xyz[0, :, :]))
        # remove deleted atoms from trajectory
        molecule = []
        for i in range(len(positions)):
            if i not in deleted_atoms:
                molecule.append(positions[i, :])
        r['molecule_xyz'] = np.array(molecule)
        r['topology'] = trunc_topology
    else:
        r['molecule_xyz'] = np.squeeze(np.array(default_traj.xyz[0, :, :]))
        r['topology'] = default_traj.topology

    r['n_residues'] = default_traj.n_residues
    table, bonds = r['topology'].to_dataframe()
    index = table.index
    r['CA'] = np.array(index[table['name'] == 'CA'])  # folding
    r['N'] = np.array(index[table['name'] == 'N'])
    r['C'] = np.array(index[table['name'] == 'C'])
    r['O'] =  np.array(index[table['name'] == 'O'])  # phi
    r['H'] =  np.array(index[table['name'] == 'H'])  # psi
    r['CB'] = []
    res_start= []
    start_resseq = table['resSeq'][0]
    end_resseq = table['resSeq'].values[-1]
    for i in range(start_resseq, end_resseq+1, 1):
        indz =  np.array(index[table['resSeq'] == i])
        res_start.append(indz[0])
    for i in range(r['n_residues']):
        if table['resName'][res_start[i]] == 'GLY':
            r['CB'].append(res_start[i]+1)
        else:
            r['CB'].append(res_start[i]+3)
    r['n_atoms'] = len(r['molecule_xyz'])
    r['noise_var_range'] = np.array([10,1])
    r['beta_cycles'] = np.array([1, cycles])
    r['beta'] = 1
    r['cycles']=cycles
    r['segment_sizes'] = np.array(segment_sizes)
    r['iterations'] = np.array(iterations)

   
    # uncomment below lines to get new data for GPU to run!
    print('n_residues: ', r['n_residues'])
    flattened_distogram = [] #np.zeros(1,r['n_residues']*r['n_residues'] * 64)
    for i in range(r['n_residues']):
        for j in range(r['n_residues']):
            flattened_distogram += list(r['distogram'][i][j][:])

    flattened_distogram = np.array(flattened_distogram)
    print('flattened distogram shape: ', flattened_distogram.shape)

    flattened_log_distogram = np.log(flattened_distogram)

    flattened_state = []
    for mol in r['molecule_xyz']:
        for axis in mol:
            flattened_state += [axis]

    print("State length: ", np.array(flattened_state).shape)


    # write to csv file, one row per variable
    stri = r'C:\Users\Anirudh Ghantasala\Documents\GitHub\Molecular_simulations\md\AlphaFold\fpga_emulator\fpga_run_data\gpudata_{}.csv'.format(name)
    f = open(stri, 'w')
    writer = csv.writer(f)
    writer.writerow(flattened_distogram)
    writer.writerow(flattened_state)
    writer.writerow(r['CA'])
    writer.writerow(r['CB'])
    writer.writerow(r['C'])
    writer.writerow(r['N'])
    # writer.writerow(r['O'])
    # writer.writerow(r['H'])
    f.close()

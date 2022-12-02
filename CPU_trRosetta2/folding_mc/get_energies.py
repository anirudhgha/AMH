# Get energies for a set of input angles for the protein T0955. Store the outputs as a pickled dataframe. 


import os,json
import numpy as np
from arguments import *
from utils import *
from pyrosetta import *
import pandas as pd
import time

vdw_weight = {0: 3.0, 1: 5.0, 2: 10.0}
rsr_dist_weight = {0: 3.0, 1: 2.0, 3: 1.0}
rsr_orient_weight = {0: 1.0, 1: 1.0, 3: 0.5}

def main():

    ########################################################
    # process inputs
    ########################################################

    # read params
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    with open(scriptdir + '/data/params_get_energies.json') as jsonfile:
        params = json.load(jsonfile)

    # get command line arguments
    args = get_args(params)
    print(args)
    # if os.path.exists(args.OUT):
    #     return

    # init PyRosetta
    init_cmd = list()
    init_cmd.append("-multithreading:interaction_graph_threads 1 -multithreading:total_threads 1")
    init_cmd.append("-hb_cen_soft")
    init_cmd.append("-detect_disulf -detect_disulf_tolerance 2.0") # detect disulfide bonds based on Cb-Cb distance (CEN mode) or SG-SG distance (FA mode)
    init_cmd.append("-relax:dualspace true -relax::minimize_bond_angles -default_max_cycles 200")
    init_cmd.append("-mute all")
    init(" ".join(init_cmd))

    # read and process restraints & sequence
    seq = read_fasta(args.FASTA)
    L = len(seq)
    params['seq'] = seq
    rst = gen_rst(params)

    ########################################################
    # Scoring functions and movers
    ########################################################
    sf = ScoreFunction()
    sf.add_weights_from_file(scriptdir + '/data/scorefxn_get_energies.wts')

    # sf1 = ScoreFunction()
    # sf1.add_weights_from_file(scriptdir + '/data/scorefxn1.wts')

    # sf_vdw = ScoreFunction()
    # sf_vdw.add_weights_from_file(scriptdir + '/data/scorefxn_vdw.wts')

    # sf_cart = ScoreFunction()
    # sf_cart.add_weights_from_file(scriptdir + '/data/scorefxn_cart.wts')

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(True)

 

    if not os.path.exists("%s/before_relax.pdb"%('.'.join(args.OUT.split('.')[:-1]))): 
        ########################################################
        # initialize pose
        ########################################################
        pose0 = pose_from_sequence(seq, 'centroid')

        # mutate GLY to ALA
        for i,a in enumerate(seq):
            if a == 'G':
                mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'ALA')
                mutator.apply(pose0)
                print('mutation: G%dA'%(i+1))

        if (args.bb == ''):
            print('not setting random (phi,psi,omega)...')
            # set_random_dihedral(pose0)
        else:
            print('setting predicted (phi,psi,omega)...')
            bb = np.load(args.bb)
            set_predicted_dihedral(pose0,bb['phi'],bb['psi'],bb['omega'])

      
      
        # The below 4 lines were originally inside the for run in range(NRUNS) loop below. 
        # They have been brough there so that sf does not change, so that alt-mh may be compared directly with lbgfs
        # sf.set_weight(rosetta.core.scoring.vdw, vdw_weight.setdefault(0, 10.0))
        sf.set_weight(rosetta.core.scoring.atom_pair_constraint, rsr_dist_weight.setdefault(0, 1.0))
        # sf.set_weight(rosetta.core.scoring.dihedral_constraint, rsr_orient_weight.setdefault(0, 0.5))
        # sf.set_weight(rosetta.core.scoring.angle_constraint, rsr_orient_weight.setdefault(0, 0.5))
            
        # 1mil samples ~ 20k seconds ~ 5.5hr
        num_samples = 100000
        num_aa = len(pose0.sequence())
        pose0.dump_pdb('dumping starting pose.pdb')

        df_angles = pd.DataFrame(columns = ['angles', 'energy'])
        df_positions = pd.DataFrame(columns = ['positions', 'Energy'])

        start = time.process_time()
        for sample in range(num_samples):
          
            pose = Pose()
            pose.assign(pose0)
            pose.remove_constraints()

            # randomize backbone torsion angles
            dphi = np.random.uniform(-180,180,num_aa)
            dpsi = np.random.uniform(-180,180,num_aa)
                        
            for i in range(1,num_aa+1):
                if pose.phi(i) + dphi[i-1] < 0:
                    print('original angle: {}, noise: {}'.format(pose.phi(i), dphi[i-1]))
                pose.set_phi(i,pose.phi(i)+dphi[i-1])
                pose.set_psi(i,pose.psi(i)+dpsi[i-1])

            add_rst(pose, rst, 3, len(seq), params)
            # pose.conformation().detect_disulfides() 
            E = sf(pose)

            reses = [i for i in range(1,num_aa+1)]

            # get pairwise distances and energy
            res1, res2, dist, e = [], [], [], []
            for i in range(1, num_aa+1):
                for j in range(1, num_aa+1):
                    res1.append(i)
                    res2.append(j)

                    # dist
                    CA1_xyz = pose.residue(i).xyz("CA")
                    CA2_xyz = pose.residue(j).xyz("CA")  
                    CA_CA_vector = CA1_xyz - CA2_xyz
                    CA_CA_vector = CA_CA_vector.norm()

                    print(i, j, CA_CA_vector)
                    # energy
            print('size', pose.energies().energies_updated())
            # a = pose.energies().residue_total_energies_array()
            # a = pose.energies().residue_pair_energies_array(reses, reses)
            summ=0
            print(len(a))
            for i in a:
                print(len(i))
                for j in i:
                    summ+=j
            print(summ)
            print(a)
            df_energies = df_energies.append({'angles':angles, 'energy': E}, ignore_index=True)

            # pose.dump_pdb('after all the changes.pdb')

            if sample % int(num_samples/100) == 0: 
                print('{}%...{} s'.format(int(sample/(num_samples)*100), time.process_time()-start))
    # if you want to read the energies
    # df_angles.to_csv('df_angles.csv')
    
    # more efficient storage
    df_angles.to_pickle('df_angles.pkl')

if __name__ == '__main__':
    main()



# THIS CODE IS LARGELY AS IS IN THE TRROSETTA2 GITHUB
# modifications are made by Lakshmi A. Ghantasala, Purdue University. 

import sys,os,json
import tempfile
import numpy as np
import time

from arguments import *
from utils import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pfold import pfold
from matplotlib import pyplot as plt
import pandas as pd

vdw_weight = {0: 3.0, 1: 5.0, 2: 10.0}
rsr_dist_weight = {0: 3.0, 1: 2.0, 3: 1.0}
rsr_orient_weight = {0: 1.0, 1: 1.0, 3: 0.5}

def main():

    ########################################################
    # process inputs
    ########################################################

    # read params
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    with open(scriptdir + '/data/params.json') as jsonfile:
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
    sf.add_weights_from_file(scriptdir + '/data/scorefxn.wts')

    sf1 = ScoreFunction()
    sf1.add_weights_from_file(scriptdir + '/data/scorefxn1.wts')

    sf_vdw = ScoreFunction()
    sf_vdw.add_weights_from_file(scriptdir + '/data/scorefxn_vdw.wts')

    sf_cart = ScoreFunction()
    sf_cart.add_weights_from_file(scriptdir + '/data/scorefxn_cart.wts')

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(True)

    min_mover1 = MinMover(mmap, sf1, 'lbfgs_armijo_nonmonotone', 0.001, True)
    min_mover1.max_iter(1000)

    min_mover_vdw = MinMover(mmap, sf_vdw, 'lbfgs_armijo_nonmonotone', 0.001, True)
    min_mover_vdw.max_iter(500)

    min_mover_cart = MinMover(mmap, sf_cart, 'lbfgs_armijo_nonmonotone', 0.000001, True)
    min_mover_cart.max_iter(300)
    min_mover_cart.cartesian(True)


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

        remove_clash(sf_vdw, min_mover_vdw, pose0)

        Emin = 99999

        ########################################################
        # minimization
        ########################################################
        score_history = []
        ref_pdb = args.REF_FOR_RMSD
        ref_pose = pyrosetta.pose_from_pdb(ref_pdb)

        # The below 4 lines were originally inside the for run in range(NRUNS) loop below. 
        # They have been brough there so that sf does not change, so that alt-mh may be compared directly with lbgfs
        sf.set_weight(rosetta.core.scoring.vdw, vdw_weight.setdefault(0, 10.0))
        sf.set_weight(rosetta.core.scoring.atom_pair_constraint, rsr_dist_weight.setdefault(0, 1.0))
        sf.set_weight(rosetta.core.scoring.dihedral_constraint, rsr_orient_weight.setdefault(0, 0.5))
        sf.set_weight(rosetta.core.scoring.angle_constraint, rsr_orient_weight.setdefault(0, 0.5))
            
        overall_min_energy, min_counter = 100000, 1
        perturb_variance_max, perturb_variance_min = 10, 1
        target_energy_reached = False
        samples = 0

        Emin = sf(pose0)
        
        pose = Pose()
        pose.assign(pose0)
        pose.remove_constraints()
        add_rst(pose, rst, 3, len(seq), params)

        for run in range(params['NRUNS']):
            # define repeat_mover here!! (update vdw weights: weak (1.0) -> strong (10.0)
            # sf.set_weight(rosetta.core.scoring.vdw, vdw_weight.setdefault(run, 10.0))
            # sf.set_weight(rosetta.core.scoring.atom_pair_constraint, rsr_dist_weight.setdefault(run, 1.0))
            # sf.set_weight(rosetta.core.scoring.dihedral_constraint, rsr_orient_weight.setdefault(run, 0.5))
            # sf.set_weight(rosetta.core.scoring.angle_constraint, rsr_orient_weight.setdefault(run, 0.5))
            
            min_mover = MinMover(mmap, sf, 'lbfgs_armijo_nonmonotone', 0.001, True)
            # min_mover_custom = MinMover()
            # min_mover_custom.score_function(sf)
            # min_mover_custom.set_movemap(mmap)
            # min_mover_custom.set_type('lbfgs_armijo_nonmonotone')
            # min_mover_custom.max_iter(100)

            # print('MIN MOVER :---------------------------------------', min_mover)
            # print('CUSTOM MIN MOVER:--------------------------', min_mover_custom)

            min_mover.max_iter(1000)

            repeat_mover = RepeatMover(min_mover, 3)

            # CHANGE BY ANIRUDH: 
            # uncomment the following 3 lines if you want pose to reset to pose0 before each noisy restart. 
            # If you want it to continue without resetting, leave it commented. 
            pose = Pose()
            pose.assign(pose0)
            pose.remove_constraints()
            add_rst(pose, rst, 3, len(seq), params)

            # noisy restarts
            # pose.assign(pose0)

            if run > 0 and (args.mode == 2 or args.mode == 6 or args.mode == 7):
                E = sf(pose)
                print(f'Energy pre-noise: {E}')
                # add_rst(pose, rst, 3, len(seq), params)

                # diversify backbone
                dphi = np.random.uniform(-10,10,L)
                dpsi = np.random.uniform(-10,10,L)
                for i in range(1,L+1):
                    pose.set_phi(i,pose.phi(i)+dphi[i-1])
                    pose.set_psi(i,pose.psi(i)+dpsi[i-1])

                # remove clashes
                # remove_clash(sf_vdw, min_mover_vdw, pose)
            else:
                E = sf(pose)
                print(f'Energy pre-noise: {E}')
            E = sf(pose)
            print(f'Energy post-noise: {E}')

            if args.mode == 0:

                # short
                print('short')
                add_rst(pose, rst, 3, 12, params)
                repeat_mover.apply(pose)
                remove_clash(sf_vdw, min_mover1, pose)
                min_mover_cart.apply(pose)
                if args.save_chk:
                    pose.dump_pdb("%s_run%d_mode%d_step%d.pdb"%('.'.join(args.OUT.split('.')[:-1]), run, args.mode, 0))

                # medium
                print('medium')
                add_rst(pose, rst, 12, 24, params)
                repeat_mover.apply(pose)
                remove_clash(sf_vdw, min_mover1, pose)
                min_mover_cart.apply(pose)
                if args.save_chk:
                    pose.dump_pdb("%s_run%d_mode%d_step%d.pdb"%('.'.join(args.OUT.split('.')[:-1]), run, args.mode, 1))

                # long
                print('long')
                add_rst(pose, rst, 24, len(seq), params)
                repeat_mover.apply(pose)
                remove_clash(sf_vdw, min_mover1, pose)
                min_mover_cart.apply(pose)
                if args.save_chk:
                    pose.dump_pdb("%s_run%d_mode%d_step%d.pdb"%('.'.join(args.OUT.split('.')[:-1]), run, args.mode, 2))

            elif args.mode == 1:

                # short + medium
                print('short + medium')
                add_rst(pose, rst, 3, 24, params)
                repeat_mover.apply(pose)
                remove_clash(sf_vdw, min_mover1, pose)
                min_mover_cart.apply(pose)
                if args.save_chk:
                    pose.dump_pdb("%s_run%d_mode%d_step%d.pdb"%('.'.join(args.OUT.split('.')[:-1]), run, args.mode, 0))

                # long
                print('long')
                add_rst(pose, rst, 24, len(seq), params)
                repeat_mover.apply(pose)
                remove_clash(sf_vdw, min_mover1, pose)
                min_mover_cart.apply(pose)
                if args.save_chk:
                    pose.dump_pdb("%s_run%d_mode%d_step%d.pdb"%('.'.join(args.OUT.split('.')[:-1]), run, args.mode, 1))
            
            # Gradient descent (lbfgs w/ noisy restarts)
            elif args.mode == 2:
                print('starting energy!!!!!!!!!!!!!!!!!', sf(pose))

                # short + medium + long
                print('short + medium + long')
                # add_rst(pose, rst, 3, len(seq), params)
                gd_iterations, minmover_iterations = 1, 100
                min_mover.max_iter(minmover_iterations)
                # print(f'Minimizing score function with Gradient Descent {minmover_iterations} times, each time for {gd_iterations} iterations')
                # pose.conformation().detect_disulfides()
                E = sf(pose)
                score_history.append(E)
                start_time = time.process_time()
                for _ in range(gd_iterations):
                    min_mover.apply(pose)
                    samples += minmover_iterations
                    E = sf(pose)
                    for gdite in range(minmover_iterations):
                        score_history.append(E)
                    if E < overall_min_energy:
                        overall_min_energy = E
                        pose.dump_pdb(args.OUT+'/pre_relaxation_gd.pdb')

                        # if E < args.T_ENERGY:
                        #     print('ener: {}, target energy: {}'.format(E, args.T_ENERGY))
                        #     target_energy_reached = True
                        #     break
                print('[gd] run: {}/{}, energy: {}, min_energy: {}, time: {}s'.format(run, params['NRUNS'], score_history[-1], overall_min_energy, time.process_time()-start_time))

                # remove_clash(sf_vdw, min_mover1, pose)
                # min_mover_cart.apply(pose)
                
                if args.save_chk:
                    pose.dump_pdb("%s_run%d_mode%d_step%d.pdb"%('.'.join(args.OUT.split('.')[:-1]), run, args.mode, 0))

            # Alternating MH
            elif args.mode == 3:
            ##############################################################################################################
            # Ani: Alternating MH 
            ##############################################################################################################
                # set these alt-mh parameters
                cycles=1
                segment_sizes = [10, 1000] 
                iterations = [20, 80]
                

                # add_rst(pose, rst, 3, len(seq), params)
                # pose.conformation().detect_disulfides()
                start_time = time.time()
                
                pfolder = pfold(pose, sf, score_history=score_history)
                pfolder.set_perturbation_variance(start=10/min_counter, end=1/min_counter)
                pfolder.set_reference_pdb(ref_pdb)
                pfolder.run_alternating_mh(cycles=cycles, iterations=iterations, segment_sizes=segment_sizes, run=run, total_runs = params['NRUNS'], beta=1)

                pose.assign(pfolder.minimum_score_pose)
                
                score_history = pfolder.score_history.copy()
                print('run {} took {} seconds'.format(run, time.time()-start_time))
                
                if pfolder.minimum_score < overall_min_energy:
                    pfolder.minimum_score_pose.dump_pdb(args.OUT+'/pre_relaxation_altmh.pdb')
                    overall_min_energy = pfolder.minimum_score
                    min_counter = 1
                else:
                    min_counter += 1    # if we can't beat a minimum, keep reducing the variance so that smaller adjustments may lead to better results. 
                # remove_clash(sf_vdw, min_mover1, pose)
                # min_mover_cart.apply(pose)
                
            # Gradient descent with localized MH
            elif args.mode == 4:
                add_rst(pose, rst, 3, len(seq), params)
                print('Gradient Descent with localized MH')
                # pose.conformation().detect_disulfides()
                print('starting energy {} for nrestart: {}'.format(sf(pose), run))

                cycles=1
                altmh_segment_sizes = [10] 
                altmh_iterations = [10]

                pfolder = pfold(pose, sf, score_history)
                pfolder.set_perturbation_variance(start=20, end=1)

                # intelligent restart
                pfolder.run_alternating_mh(run=run, total_runs=params['NRUNS'], cycles=cycles, iterations=altmh_iterations, segment_sizes=altmh_segment_sizes, beta=1)

                pose.assign(pfolder.minimum_score_pose) #WAS USED FOR PAPER
                # pose.assign(pfolder.pose)

                # optimization
                min_mover_iterations = 90
                min_mover.max_iter(min_mover_iterations)
                min_mover.apply(pose)
                
                samples += min_mover_iterations + sum(altmh_iterations)
                E = sf(pose)
                print('ener: {}, target energy: {}'.format(E, args.T_ENERGY))
                score_history = pfolder.score_history.copy()
                for i in range(min_mover_iterations):
                    score_history.append(E)
                if E < overall_min_energy:
                    overall_min_energy = E
                
            # MH with annealing
            elif args.mode == 5:
                
                # set these mh parameters
                iterations = 100
                samples = iterations
                beta_start, beta_end = 0.1, 1
                beta = (beta_end-beta_start)/params['NRUNS'] * run #beta keeps increasing each noisy restart (i.e temp decreases)

                # add_rst(pose, rst, 3, len(seq), params)
                # pose.conformation().detect_disulfides()
                start_time = time.time()
                pfolder = pfold(pose, sf, score_history)
                pfolder.set_perturbation_variance(start=10/min_counter, end=1/min_counter)
                pfolder.set_reference_pdb(ref_pdb)
                pfolder.run_sf_mh(run=run, total_runs=params['NRUNS'], iterations=iterations, beta=beta, printlog=True)
                pfolder.save_outputs(folder=args.OUT)
                pose.assign(pfolder.minimum_score_pose)

                score_history = pfolder.score_history.copy()
                print('run {} took {} seconds'.format(run, time.time()-start_time))
                
                if pfolder.minimum_score < overall_min_energy:
                    pfolder.minimum_score_pose.dump_pdb(args.OUT+'/pre_relaxation_altmh.pdb')
                    overall_min_energy = pfolder.minimum_score
                    min_counter = 1
                else:
                    min_counter += 1    # if we can't beat a minimum, keep reducing the variance so that smaller adjustments may lead to better results. 
                # remove_clash(sf_vdw, min_mover1, pose)
                
               
            # MALA with L-bfgs w/ noisy restarts
            elif args.mode == 6:
                print('L-bfgs with MH criterion acceptance')
                # pose.conformation().detect_disulfides()
                print('starting energy {} for nrestart: {}'.format(sf(pose), run))
                # add_rst(pose, rst, 3, len(seq), params)

                # Generate gradient-based proposal
                pose_proposal= Pose()
                add_rst(pose_proposal, rst, 3, len(seq), params)


                min_pose = Pose()
                min_pose.assign(pose)
                min_mover_iterations = 20
                min_mover.max_iter(min_mover_iterations)
                for i in range(100):
                    pose_proposal.assign(pose)

                    min_mover.apply(pose_proposal)
                    samples += min_mover_iterations
                     
                    di = sf(pose)
                    df = sf(pose_proposal)
                    Pf = np.exp(-1 * df)
                    Pi = np.exp(-1 * di)
                    f = np.abs(Pf / Pi) - np.random.rand()
                        
                    # accept or reject
                    if df < di or f > 0:
                        pose.assign(pose_proposal) 
                        di = df

                    score_history.append(di)
                    if di < overall_min_energy:
                        # print("New Min found! (iter=%d): %.1f --> %.1f"%(run, Emin, di))
                        overall_min_energy = di
                        min_pose.assign(pose)
                pose.assign(min_pose)
            
            # single flip mh w/ noisy restarts
            elif args.mode == 7:
                
                # set these mh parameters
                iterations = 100
                samples = iterations
                beta_start, beta_end = 0.1, 1
                beta = (beta_end-beta_start)/params['NRUNS'] * run #beta keeps increasing each noisy restart (i.e temp decreases)

                # add_rst(pose, rst, 3, len(seq), params)
                # pose.conformation().detect_disulfides()
                start_time = time.time()
                pfolder = pfold(pose, sf, score_history)
                pfolder.set_perturbation_variance(start=10/min_counter, end=1/min_counter)
                pfolder.set_reference_pdb(ref_pdb)
                pfolder.run_sf_mh(run=run, total_runs=params['NRUNS'], iterations=iterations, beta=beta, printlog=True)
                pfolder.save_outputs(folder=args.OUT)
                pose.assign(pfolder.minimum_score_pose)

                score_history = pfolder.score_history.copy()
                print('run {} took {} seconds'.format(run, time.time()-start_time))
                
                if pfolder.minimum_score < overall_min_energy:
                    pfolder.minimum_score_pose.dump_pdb(args.OUT+'/pre_relaxation_altmh.pdb')
                    overall_min_energy = pfolder.minimum_score
                    min_counter = 1
                else:
                    min_counter += 1    # if we can't beat a minimum, keep reducing the variance so that smaller adjustments may lead to better results. 
                # remove_clash(sf_vdw, min_mover1, pose)
                
               
            # check whether energy has decreased
            # pose.conformation().detect_disulfides() # detect disulfide bonds
            
            E = sf(pose)
            print(f'energy after 100 iterations: {E}, current min: {Emin}')
            # print("Energy(iter=%d): %.1f (accept)"%(run, E))
            # pose0.assign(pose)

            if E < Emin:
                print("New Min found! (iter=%d): %.1f --> %.1f (accept)"%(run, Emin, E))
                Emin = E
                pose0.assign(pose)
            # else:
                # print("Energy(iter=%d): %.1f --> %.1f (reject)"%(run, Emin, E))
            rmsd_to_ref = pyrosetta.rosetta.core.scoring.CA_rmsd(pose0, ref_pose)
            print('minimum energy (iter:{}) RMSD to reference {}: {}'.format(run, ref_pdb, rmsd_to_ref))



        # mutate ALA back to GLY
        for i,a in enumerate(seq):
            if a == 'G':
                mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'GLY')
                mutator.apply(pose0)
                print('mutation: A%dG'%(i+1))
        
        # dump information after minimization
        if Emin == 99999:
            Emin = sf(pose)
        
        info = {'mode': args.mode, 'min_energy': Emin, 'samples': samples, 'rmsd_to_ref':rmsd_to_ref, 'target_energy': args.T_ENERGY, 'score_history':score_history}
        df = pd.DataFrame(info)
        df.to_csv(args.OUT+'/info.csv')
        pose0.dump_pdb(args.OUT+'/pre_relaxation.pdb')
        plt.plot(score_history)
        plt.title('Energy during Protein Optimization')
        plt.savefig(args.OUT+'/energy_plot.pdf')
            



        ########################################################
        # fix backbone geometry
        ########################################################
        pose0.remove_constraints()

        # apply more strict criteria to detect disulfide bond
        # Set options for disulfide tolerance -> 1.0A
        print (rosetta.basic.options.get_real_option('in:detect_disulf_tolerance'))
        rosetta.basic.options.set_real_option('in:detect_disulf_tolerance', 1.0)
        print (rosetta.basic.options.get_real_option('in:detect_disulf_tolerance'))
        pose0.conformation().detect_disulfides()

        # Converto to all atom representation
        switch = SwitchResidueTypeSetMover("fa_standard")
        switch.apply(pose0)

        # idealize problematic local regions if exists
        idealize = rosetta.protocols.idealize.IdealizeMover()
        poslist = rosetta.utility.vector1_unsigned_long()

        scorefxn=create_score_function('empty')
        scorefxn.set_weight(rosetta.core.scoring.cart_bonded, 1.0)
        scorefxn.score(pose0)

        emap = pose0.energies()
        print("idealize...")
        for res in range(1,L+1):
            cart = emap.residue_total_energy(res)
            if cart > 50:
                poslist.append(res)
                print( "idealize %d %8.3f"%(res,cart) )

        if len(poslist) > 0:
            idealize.set_pos_list(poslist)
            try:
                idealize.apply(pose0)

            except:
                print('!!! idealization failed !!!')

        # Save checkpoint
        if args.save_chk:
            pose0.dump_pdb("%s_before_relax.pdb"%'.'.join(args.OUT.split('.')[:-1]))

    else: # checkpoint exists
        pose0 = pose_from_file("%s/before_relax.pdb"%('.'.join(args.OUT.split('.')[:-1])))
        #
        print ("to centroid")
        switch_cen = SwitchResidueTypeSetMover("centroid")
        switch_cen.apply(pose0)
        #
        print ("detect disulfide bonds")
        # Set options for disulfide tolerance -> 1.0A
        print (rosetta.basic.options.get_real_option('in:detect_disulf_tolerance'))
        rosetta.basic.options.set_real_option('in:detect_disulf_tolerance', 1.0)
        print (rosetta.basic.options.get_real_option('in:detect_disulf_tolerance'))
        pose0.conformation().detect_disulfides()
        #
        print ("to all atom")
        switch = SwitchResidueTypeSetMover("fa_standard")
        switch.apply(pose0)


    ########################################################
    # full-atom refinement
    ########################################################

    if args.fastrelax == True:
        mmap = MoveMap()
        mmap.set_bb(True)
        mmap.set_chi(True)
        mmap.set_jump(True)
        
        # First round: Repeat 2 torsion space relax w/ strong disto/anglogram constraints
        sf_fa_round1 = create_score_function('ref2015_cart')
        sf_fa_round1.set_weight(rosetta.core.scoring.atom_pair_constraint, 3.0)
        sf_fa_round1.set_weight(rosetta.core.scoring.dihedral_constraint, 1.0)
        sf_fa_round1.set_weight(rosetta.core.scoring.angle_constraint, 1.0)
        sf_fa_round1.set_weight(rosetta.core.scoring.pro_close, 0.0)
        
        relax_round1 = rosetta.protocols.relax.FastRelax(sf_fa_round1, "%s/data/relax_round1.txt"%scriptdir)
        relax_round1.set_movemap(mmap)
        
        print('relax: First round... (focused on torsion space relaxation)')
        params['PCUT'] = 0.15
        pose0.remove_constraints()
        add_rst(pose0, rst, 3, len(seq), params, nogly=True, use_orient=True)
        relax_round1.apply(pose0)
       
        # Set options for disulfide tolerance -> 0.5A
        print (rosetta.basic.options.get_real_option('in:detect_disulf_tolerance'))
        rosetta.basic.options.set_real_option('in:detect_disulf_tolerance', 0.5)
        print (rosetta.basic.options.get_real_option('in:detect_disulf_tolerance'))

        sf_fa = create_score_function('ref2015_cart')
        sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 0.1)
        sf_fa.set_weight(rosetta.core.scoring.dihedral_constraint, 0.0)
        sf_fa.set_weight(rosetta.core.scoring.angle_constraint, 0.0)
        
        relax_round2 = rosetta.protocols.relax.FastRelax(sf_fa, "%s/data/relax_round2.txt"%scriptdir)
        relax_round2.set_movemap(mmap)
        relax_round2.cartesian(True)
        relax_round2.dualspace(True)

        print('relax: Second round... (cartesian space)')
        params['PCUT'] = 0.30 # To reduce the number of pair restraints..
        pose0.remove_constraints()
        pose0.conformation().detect_disulfides() # detect disulfide bond again w/ stricter cutoffs
        # To reduce the number of constraints, only pair distances are considered w/ higher prob cutoffs
        add_rst(pose0, rst, 3, len(seq), params, nogly=True, use_orient=False, p12_cut=params['PCUT'])
        # Instead, apply CA coordinate constraints to prevent drifting away too much (focus on local refinement?)
        add_crd_rst(pose0, L, std=1.0, tol=2.0)
        relax_round2.apply(pose0)

        # Re-evaluate score w/o any constraints
        scorefxn_min=create_score_function('ref2015_cart')
        scorefxn_min.score(pose0)

    ########################################################
    # save final model
    ########################################################
    pose0.dump_pdb(args.OUT+'/final.pdb')
    # print rmsd to reference structure
    ref_pose = pyrosetta.pose_from_pdb(ref_pdb)
    rmsd_to_ref = pyrosetta.rosetta.core.scoring.CA_rmsd(pose0, ref_pose)
    print('final RMSD to reference {}: {}'.format(ref_pdb, rmsd_to_ref))


if __name__ == '__main__':
    main()

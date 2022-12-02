import pyrosetta
from pyrosetta import *
import numpy as np
import matplotlib.pyplot as plt
from random import gauss
import time
from utils import *



# from numba import jit
# import mdtraj as md
# import csv
# import MDAnalysis.analysis.rms as rms

    # [To Add] Read in rmsd of trRosetta code vs native protein
    # [To Add] Read in energy for best solution of trRosetta default code to compare against
    # [To Add] Read in reference pose for real protein to calculate rmsd against
    # [To Add] save_outputs: csv information
    # [To Add] save_outputs: RMSD vs proposals plot
    # [To Add] save a trajectory every x number of frames, other metrics

class pfold:
    def __init__(self, pose, scorefxn, score_history = []):
        self.pose = pose
        self.scorefxn = scorefxn
        self.n_residues = pose.size()

        self.segment_size, self.iterations = 0, 0
        self.noise_var_range = [0,0]
        self.num_proposals = 0
        self.accept_count = 0
        self.minimum_score_pose = Pose()
        self.minimum_score_pose.assign(self.pose)# = self.pose.clone()
        self.minimum_score= self.scorefxn(self.pose)
        self.score_history=score_history[:]
        self.track_min_energy = []
        self.reference_pdb = ''
        self.beta = 1 
        self.annealing_beta_start = 1
        self.annealing_beta_end = 1
        self.proposals_to_min_score = 0
        self.annealing = False

    def reinitialize_pose(self, pose):
        self.pose = pose
        
    def set_annealing(self, beta_start, beta_end):
        self.annealing_beta_start = beta_start
        self.annealing_beta_end = beta_end
        self.annealing = True

    def set_beta(self, beta):
        self.beta = beta

    def set_perturbation_variance(self, start, end):
        self.noise_var_range = np.array([start, end])

    def set_reference_pdb(self, reference_pdb):
        self.reference_pdb = reference_pdb
        self.ref_pose = pyrosetta.pose_from_pdb(self.reference_pdb)

    def save_outputs(self, folder, print_outputs=True):
        self.acceptance_prob = self.accept_count / self.num_proposals
        self.rmsd = pyrosetta.rosetta.core.scoring.CA_rmsd(self.pose, self.ref_pose)
        # write final structure pdb
        # self.minimum_score_pose.dump_pdb("{}/min_energy_altmh.pdb".format(folder))
        # Main plot (Energy vs proposals with minimum energy rmsd tracking)
        # plt.figure(0)
        # plt.title('Energy during M-H Iterations on Protein Chain')
        # plt.plot(self.score_history, color='black', label='energy')
        # plt.legend()
        # plt.savefig('{}/energy_plot.pdf'.format(folder))
        # save score history
        # np.savetxt("{}/e_dump.csv".format(self.score_history), delimiter=",")
        # with open(folder+'/info.txt', 'w') as f:
        #     f.write('[Alt-mh log] Minimum Energy: {}\n'.format(self.minimum_score))
        #     f.write('[Alt-mh log] proposals to minimum score (pre-relax): {}\n'.format(self.proposals_to_min_score))
        #     f.write('[Alt-mh log] RMSD of min_score pose to reference: {}\n'.format(self.rmsd))
        #     f.write('[Alt-mh log] Proposals: {}\n'.format(self.num_proposals))
        #     f.write('[Alt-mh log] Acceptance %: {}\n'.format(self.acceptance_prob))

        if print_outputs:            
            print('[Alt-mh log] Proposals: {}'.format(self.num_proposals))
            print('[Alt-mh log] Acceptance %: {}'.format(self.acceptance_prob))
            print('[Alt-mh log] Minimum Energy: {}'.format(self.minimum_score))
            print('[Alt-mh log] RMSD of min_score pose to reference: ', self.rmsd)
            print('[Alt-mh log] proposals to minimum score (pre-relax): {}\n'.format(self.proposals_to_min_score))


    def var_bound(val, lower, upper):
        if val < lower:
            return lower
        elif val > upper:
            return upper
        else:
            return val
    
    def perturb_torsion(pose, residue, angle, variance):
        delta = gauss(0, variance)
        if angle is 'phi':
            pose.set_phi(residue, pose.phi(residue) + delta)
        else:
            pose.set_psi(residue, pose.psi(residue) + delta)
        return pose

    def get_metric_energy(self, pose):
        return self.metric_scorefxn(pose)
    
    # Run alternating metropolis-hastings to minimize scorefxn. Assumes that restraints are already added to self.pose
    def run_alternating_mh(self, cycles, iterations, segment_sizes, run=1, total_runs=0, beta=None, printlog=True):

        proposal_pose = Pose()
        # add_rst(proposal_pose, rst, 3, len(seq), params)

        di, df, Pi, Pf = 0, 0, 0, 0
        di = self.scorefxn(self.pose)
        print('initial energy: ', di)

        residue_mask_initialize = [False for _ in range(self.n_residues)]
        self.minimum_score = self.scorefxn(self.pose)
        noise_var_delta = (self.noise_var_range[0] - self.noise_var_range[1]) / (cycles * np.sum(iterations) * (self.n_residues) * 2)
        noise_var = self.noise_var_range[0]
        total_iterations = cycles*sum(iterations) * total_runs
        iteration_counter = run*cycles*sum(iterations)
        
        

        if beta:
            self.beta = beta

        # precompute the residue masks, one mask per segment per residue (use like mask_dict[segment][res])
        mask_dict = {}
        for select in range(len(iterations)):
            segment_size = segment_sizes[select]
            masks_res_dict = {}
            for res in range(1, self.n_residues+1):

                start = pfold.var_bound(res - int(segment_size/2), 1, self.n_residues)
                end = pfold.var_bound(res + int(segment_size/2), 1, self.n_residues)

                mask = residue_mask_initialize[:]
                for m in range(start-1, end):
                    mask[m] = True
                residue_mask = pyrosetta.rosetta.utility.vector1_bool() 
                residue_mask.extend(mask)

                masks_res_dict[res] = residue_mask
            mask_dict[segment_size] = masks_res_dict

        for cycle in range(cycles):
            for select in range(len(iterations)):
                self.iterations = iterations[select]
                self.segment_size = segment_sizes[select]

                for ite in range(self.iterations):
                    prev_accept_count = self.accept_count

                    # precompute perturbations to save time 
                    noise = np.random.normal(0, noise_var, (self.n_residues, 2))

                    for res in range(1, self.n_residues+1):
                       
                        # start and end for localized energy calculation, start, end go from 1-->n_residue, not 0-->n_residue-1
                        start = pfold.var_bound(res - int(self.segment_size/2), 1, self.n_residues)
                        end = pfold.var_bound(res + int(self.segment_size/2), 1, self.n_residues)

                        for angle in ['phi', 'psi']:

                            # generate proposal
                            proposal_pose.assign(self.pose)
                            if angle is 'phi':
                                delta = noise[res-1,0]
                                proposal_pose.set_phi(res, proposal_pose.phi(res) + delta)
                            else:
                                delta = noise[res-1,1]
                                proposal_pose.set_psi(res, proposal_pose.psi(res) + delta)

                            # score proposal and pose
                            self.scorefxn(proposal_pose) 
                            df = self.scorefxn.get_sub_score(proposal_pose, mask_dict[self.segment_size][res])
                            # print('proposal score: ', df)
                            self.scorefxn(self.pose)
                            di = self.scorefxn.get_sub_score(self.pose, mask_dict[self.segment_size][res]) 
                            Pf = np.exp(-1 * df * beta)
                            Pi = np.exp(-1 * di * beta)

                            # MH critereon
                            f = np.abs(Pf / Pi) - np.random.rand()
                            
                            # accept and track metrics
                            if df < di or f > 0:
                                self.pose.assign(proposal_pose) 
                                self.accept_count += 1
                                full_score = proposal_pose.energies().total_energy()

                                if full_score < self.minimum_score:
                                    self.minimum_score = full_score
                                    self.minimum_score_pose.assign(self.pose) # = self.pose.clone()
                                    self.proposals_to_min_score = self.num_proposals
                
                            # anneal perturbation variance
                            noise_var -= noise_var_delta

                            self.num_proposals += 1

                    iteration_counter += 1
                    self.score_history.append(full_score)
                    print('[pfold] ite: {}/{}, segment: {}, accepted: {}, noise_var: {}, beta: {}, min_energy: {}, current_energy: {}'.format(iteration_counter, total_iterations, self.segment_size, self.accept_count - prev_accept_count, round(noise_var,2), round(self.beta, 2), round(self.minimum_score,2), round(self.score_history[-1],2)))


    def run_sf_mh(self, iterations, beta=None, total_runs=0, run=1, printlog=True):

        proposal_pose = Pose()
        df, Pi, Pf = 0, 0, 0
        di = self.scorefxn(self.pose)
        self.minimum_score = self.scorefxn(self.pose)
        self.iterations = iterations
        total_iterations = total_runs * iterations
        self.score_history.append(di)

        if beta:
            self.beta = beta

        noise_var_delta = (self.noise_var_range[0] - self.noise_var_range[1]) / (iterations * self.n_residues * 2)
        noise_var = self.noise_var_range[0]

        

        if printlog:
            print('[pfold] Starting sf_mh: iterations: {}, perturbation_var: {}, starting energy: {}'.format(self.iterations, round(noise_var, 2), self.minimum_score))

        for ite in range(self.iterations):
            prev_accept_count = self.accept_count

            # precompute perturbations to save time 
            noise = np.random.normal(0, noise_var, (self.n_residues, 2))

            for res in range(1, self.n_residues+1):
                
                for angle in ['phi', 'psi']:

                    # gen proposal
                    start = time.process_time()
                    proposal_pose.assign(self.pose)
                    delta = gauss(0, noise_var)
                    if angle is 'phi':
                        delta = noise[res-1,0]
                        proposal_pose.set_phi(res, proposal_pose.phi(res) + delta)
                    else:
                        delta = noise[res-1,1]
                        proposal_pose.set_psi(res, proposal_pose.psi(res) + delta)                    
                    # print('\t[time] generating proposal: {}'.format(time.process_time()-start))

                    # score proposal
                    start = time.process_time()
                    df = self.scorefxn(proposal_pose)
                    # print('\t[time] scoring proposal: {}'.format(time.process_time()-start))

                    Pf = np.exp(-1 * df * self.beta)
                    Pi = np.exp(-1 * di * self.beta)

                    # print('res {} {}: df:{} di:{}'.format(res, angle, df, di))
                    # print('pf: {} pi: {}'.format(Pf, Pi))

                    # accept proposal if f > 0
                    f = np.abs(Pf / Pi) - np.random.rand()
                    
                    # accept and track metrics
                    start = time.process_time()
                    if df < di or f > 0:
                        self.pose.assign(proposal_pose) 
                        di = df
                        self.accept_count += 1
                        full_score = proposal_pose.energies().total_energy()

                        if full_score < self.minimum_score:
                            self.minimum_score = full_score
                            self.minimum_score_pose = proposal_pose.clone()
                            self.proposals_to_min_score = self.num_proposals
                    else:
                        full_score = self.pose.energies().total_energy()
                    # print('\t[time] acceptance and tracking metrics: {}'.format(time.process_time()-start))

                    # anneal perturbation variance
                    noise_var -= noise_var_delta

                    self.num_proposals += 1

            # once per iteration        
            self.score_history.append(full_score)
            print('[pfold] ite: {}/{}, accepted: {}, noise_var: {}, beta: {}, min_energy: {}, current_energy: {}'.format(ite+(run*iterations),total_iterations, self.accept_count - prev_accept_count, round(noise_var,2), round(self.beta, 2), round(self.minimum_score,2), round(self.score_history[-1],2)))
                    
        print('[pfold] sf_mh complete')




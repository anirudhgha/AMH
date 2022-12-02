import pyrosetta
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta import pose_from_sequence, Pose
from pyrosetta.rosetta.protocols.loops import get_cen_scorefxn, get_fa_scorefxn
from pyrosetta import *
from random import random,randint,gauss


import sys,os,json
import tempfile
import numpy as np

from arguments import *
from utils import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover


pyrosetta.init()

# show results on pymol server
pmm = PyMOLMover()
pmm.keep_history(True)


seq = 'SQETRKKCTEMKKKFKNCEVRCDESNHCVEVRCSDTKYTLC'
pose0 = pose_from_sequence(seq, 'centroid')

scorefxn = get_cen_scorefxn()
# pose_from_pdb(p, “mc_initial.pdb”)
# #set up score function
# scorefxn = ScoreFunction()
# scorefxn.set_weight(hbond_sr_bb,1.0)
# scorefxn.set_weight(vdw, 1.0)
#set up MonteCarlo object
mc = MonteCarlo(pose0, scorefxn, 1.0)
#set up mover
def perturb_bb(pose, var=10):
    # propose = Pose()
    # propPose.assign(pose)
    res = randint(1, pose.size())
    delta = gauss(0, var)
    if random() < 0.5:
        pose.set_phi(res, pose.phi(res)+delta)
    else:
        pose.set_psi(res, pose.psi(res)+delta)
    return pose
#set up protocol
def my_protocol(mc, pose):
    mc.reset(pose)
    for i in range(1,1000):
        perturb_bb(pose)
        mc.boltzmann(pose)

        pmm.apply(pose)
        if (i%100 == 0):
            mc.recover_low(pose)
    #output lowest-energy structure
    mc.recover_low(pose)
    return pose

my_protocol(mc, pose0)
pose0.dump_pdb('mc_final.pdb') 

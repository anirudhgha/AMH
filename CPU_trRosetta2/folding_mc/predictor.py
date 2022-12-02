from pyrosetta.rosetta import *
from pyrosetta.toolbox import *
from math import exp
from random import random,randint,gauss
from pyrosetta import init, pose_from_pdb, SwitchResidueTypeSetMover
from pyrosetta.rosetta.protocols.loops import get_cen_scorefxn, get_fa_scorefxn
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import ClassicFragmentMover

init()

class Predictor(object):
  def __init__(self,
               applyPyMOL = False,
               nIters = 1000,
               minimize = False,
               minFreq = 10,
               name = "",
               full_atom = True,
               ):
    self.forcefield = get_fa_scorefxn() if full_atom else get_cen_scorefxn()
    if applyPyMOL:
      self.pmm = PyMOLMover()
      self.pmm.keep_history(True)
    else:
      self.pmm = None
    if minimize:
      self.minm = MinMover()
      movemap = MoveMap()
      movemap.set_bb(True)
      self.minm.movemap(movemap)
      self.minm.score_function(self.forcefield)
      self.minFreq = minFreq
    else:
      self.minm = None

    self.nIters = nIters
    self.dumpFreq = self.nIters/10
    self.name = "%s-%s%s"%("fa" if full_atom else "centroid","min-" if minimize else "",name)

  def predict(self, start):
    raise NotImplementedError("Predictor should not be instantiated")

  def sampleMove(self, current, **kwargs):
    raise NotImplementedError("Predictor should not be instantiated")


class MonteCarloPredictor(Predictor):

  # Method:  Accept Move
  # --------------------
  # Return True to accept the move and False to reject
  # according to the Metropolis Criterion
  #
  # Assume kT = 1.0
  #
  def acceptMove(self, currE, newE):
    ### BEGIN YOUR CODE HERE ###

    accept = abs(exp(-1*newE) / exp(-1*currE)) - random()
    return (accept > 0)

  ###  END YOUR CODE HERE  ###

  # Method:  Standard Monte Carlo Method
  # ------------------------------------
  # Accepts a starting pose, start, and use a Monte Carlo
  # algorithm to predict the most likely pose, optPose.
  # 
  # Useful methods:
  # self.forcefield(pose) -- evaluates "energy" of pose
  # pose.assign(newPose) -- assigns newPose's coordinates on pose
  # self.sampleMove(pose) -- returns new pose after sampling a move
  # self.acceptMove(currE, newE) -- returns True if and only if move should be accepted
  #
  def predict(self, start):
    currPose = Pose()
    currPose.assign(start)
    optPose = Pose()
    optPose.assign(currPose)

    currE = self.forcefield(currPose)
    optE = currE

    # main Monte Carlo Loop
    for i in range(self.nIters):
      


      ### BEGIN YOUR CODE HERE ###

      # use self.sampleMove(), self.forcefield, self.acceptMove(), and assign() to 
      # update the variables pose, currE, newE, newPose optPose, and optE as appropriate
      # You should be assigning a newPose based off pose and determining whether or not to 
      # accept the newPose based on the pose and newPose energies. If your newPose has the 
      # lowest energy you've seen thus far, then you should update the optimum position and energy

      ## ~8 lines of code ## 
      
      propPose = self.sampleMove(currPose)
      propE = self.forcefield(propPose)
      if self.acceptMove(currE, propE):
        currE = propE
        currPose.assign(propPose)
      if currE < optE:
        optE = currE
        optPose.assign(currPose)


      ###  END YOUR CODE HERE  ###

      

      if self.pmm:
        self.pmm.apply(currPose)
      if i%self.dumpFreq == 0:
        print("%s dumping pdb at iteration %d"%(self.name,i))
        optPose.dump_pdb("out/%s-%d-%d.pdb"%(self.name,self.nIters,i))
      if self.minm and i%self.minFreq == 0:
        self.minm.apply(currPose)
    optPose.dump_pdb("out/%s-%d-final.pdb"%(self.name,self.nIters))
    return optPose


class DihedralPredictor(MonteCarloPredictor):

  def __init__(self,
               applyPyMOL = False,
               nIters = 1000,
               minimize = False,
               name = "dihedral",
               full_atom = True,
               variance = 25,
               ):
    self.variance = variance
    super(DihedralPredictor,self).__init__(applyPyMOL = applyPyMOL,
                                           nIters = nIters,
                                           minimize = minimize,
                                           name = name,
                                           full_atom = full_atom,
                                           )

  # Method:  sample dihedral move
  # -----------------------------
  # Samples one of three styles of dihedral moves uniformly
  # at random.  Either increments phi by delta, psi by delta,
  # or phi by delta and psi by -delta.
  # Returns new pose after move.
  # 
  # Useful methods:
  # pose.phi(i) -- returns phi of residue i
  # pose.psi(i) -- returns psi of residue i
  # pose.set_phi(i,newPhi) -- sets phi_i to be newPhi
  # pose.set_psi(i,newPsi) -- sets psi_i to be newPsi
  #
  def sampleMove(self, pose):
    # newPose = Pose()
    # newPose.assign(pose)
    # # select a random residue
    # res = randint(1,len(newPose.sequence()))
    # # select the size of the deviation
    # delta = gauss(0,self.variance)
    # # randomly select the style of dihedral move
    # r = random()
    # if r < 0.33333:
    #   # dihedral move in phi
    #   muPhi = newPose.phi(res)
    #   newPose.set_phi(res,muPhi + delta)

    # elif r < 0.66667:
    #   # dihedral move in psi 
    #   muPsi = newPose.psi(res)
    #   newPose.set_psi(res,muPsi + delta)

    # else:
    #   # shear move
    #   # update phi <- phi + delta and psi <- psi - delta
    #   muPhi = newPose.phi(res)
    #   muPsi = newPose.psi(res)

    #   newPose.set_phi(res, muPhi + delta)
    #   newPose.set_psi(res, muPsi - delta)

    newPose = Pose()
    newPose.assign(pose)

    # select a random dihedral residue
    res = randint(1,len(newPose.sequence()))
    delta = gauss(0,self.variance)
    # select a random dihedral residue
    r = random()
    if r < 0.5:
      # dihedral move in phi
      newPose.set_phi(res,newPose.phi(res) + delta)

    else:
      # dihedral move in psi 
      newPose.set_psi(res,newPose.psi(res) + delta)

    return newPose


class FragmentPredictor(MonteCarloPredictor):

  def __init__(self,
               applyPyMOL = False,
               nIters = 100,
               minimize = False,
               name = "frag",
               full_atom = True,
               fragSize = 9,
               fragFile = "",
               ):
    super(FragmentPredictor,self).__init__(applyPyMOL = applyPyMOL,
                                           nIters = nIters,
                                           minimize = minimize,
                                           name = "%s%d"%(name,fragSize),
                                           full_atom = full_atom,
                                           )
    fragset = ConstantLengthFragSet(fragSize)
    fragset.read_fragment_file(fragFile)
    movemap = MoveMap()
    movemap.set_bb(True)
    self.fragMover = ClassicFragmentMover(fragset,movemap)

  # Method: sample fragment move
  # ----------------------------
  # Samples a move by replacing the coordinates for a randomly
  # selected fragment of the pose with those in the fragment
  # set.
  #
  # Useful Methods:
  # fragMover.apply(pose)
  def sampleMove(self,pose):
    newPose = Pose()
    newPose.assign(pose) 
    self.fragMover.apply(newPose)
    return newPose



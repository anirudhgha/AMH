import sys,os

if len(sys.argv) < 4:
  print("Usage: python predict.py <input>.pdb <output>.pdb -predictor [additional options]")
  sys.exit(1)

if not os.path.isfile(sys.argv[1]):
  print("%s is not a valid path to a pdb file"%sys.argv[1])
  sys.exit(1)

predictors = set(("-dihedral", "-frag9", "-frag3", "-all"))
options = set(sys.argv[3:])

try:
  predictorName = predictors.intersection(options).pop()
  options.remove(predictorName)
except KeyError:
  print("Make sure you've specified a predictor out of:")
  print([p for p in predictors])
  sys.exit(1)

pymol = "-pymol" in options
if pymol:
  options.remove("-pymol")

minimize = "-min" in options
if minimize:
  options.remove("-min")

fa = "-fa" in options
if fa:
  options.remove("-fa")

var = "-var" in options
if minimize:
  options.remove("-var")

nIters = 1000
for o in options:
  try:
    nIters = int(o[1:])
    break
  except ValueError:
    continue


from predictor import *


pose = pose_from_pdb(sys.argv[1])
pose.pdb_info().name("struct")

if not fa:
  toCentroid = SwitchResidueTypeSetMover('centroid')
  toCentroid.apply(pose)

if predictorName == "-dihedral":
  P = DihedralPredictor(applyPyMOL = pymol, nIters = nIters, full_atom = fa, variance=var)
elif predictorName == "-frag9":
  P = FragmentPredictor(applyPyMOL = pymol, nIters = nIters, fragSize = 9, fragFile = "frags/frag9", full_atom = fa)
elif predictorName == "-frag3":
  P = FragmentPredictor(applyPyMOL = pymol, nIters = nIters, fragSize = 3, fragFile = "frags/frag3", full_atom = fa)

pose = P.predict(pose)
# toFA = SwitchResidueTypeSetMover('fa_standard')
# toFA.apply(pose)
pose.dump_pdb("out/%s"%sys.argv[2])




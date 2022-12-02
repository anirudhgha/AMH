import pyrosetta
from pyrosetta import *
import rosetta
import numpy as np
import csv

init()


seq = [
        ['SQETRKKCTEMKKKFKNCEVRCDESNHCVEVRCSDTKYTLC','T0955.pdb','T0955'],                                          #_T0955
        ['GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS','6f45.pdb','T0953s1'],          #_T0953s1
        ['MSYDYSSLLGKITEKCGTQYNFAIAMGLSERTVSLKLNDKVTWKDDEILKAVHVLELNPQDIPKYFFNAKVH','6tri.pdb','T0974s1'],          #_T0974s
        ['DIYGDEITAVVSKIENVKGISQLKTRHIGQKIWAELNILVDPDSTIVQGETIASRVKKALTEQIRDIERVVVHFEPARK','6qek.pdb','T1006'],     #_T1006
        ['TDELLERLRQLFEELHERGTEIVVEVHINGERDEIRVRNISKEELKKLLERIREKIEREGSSEVEVNVHSGGQTWTFNEK','6msp.pdb','T1008'],    #_T1008
        ['GSMGKDEALEKDLNDVSKEINLMLSTYAKLLSERAAVDASYIDEIDELFKEANAIENFLIQKREFLRQR','6hk8.pdb','T1072s2'],             #_T1072s2
        ['MNVDPHFDKFMESGIRHVYMLFENKSVESSEQFYSFMRTTYKNDPCSSDFECIERGAEMAQSYARIMNIKLETE','6px4.pdb','T1046s1'],        #_T1046s1 
        ['RGPSNGQSVLENSVQVKETSPRRVSVDPQTGEFVVFDRTLGDVYHGHVRAWKDLTSDMQNALVRGGYVDRKGNPK','5t87.pdb','T0884']          #_T0884
    ] 




with open('../out/min_rmsds.csv', 'w+', newline='') as csvfile:
    w = csv.writer(csvfile, delimiter=',')

    w.writerow(['protein','1','2','3','4','5','6','7','8','9','10'])
    for s in seq:
        rmsds = [s[2]]
        for run in range(10):
            print(f'[log] ######################################################PROCESSING: {s[2]}\n\n\n\n')
            mypose = pyrosetta.pose_from_pdb(f'../proteins/{s[2]}/{run}/pre_relaxation.pdb')

            # mutate ALA back to GLY
            for i,a in enumerate(s[0]):
                if a == 'G':
                    mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'GLY')
                    mutator.apply(mypose)
                    print('mutation: A%dG'%(i+1))    
                # pose.dump_pdb(f'../out/{protein}/mutated_back_alt.pdb')

            # find rmsd between pose and original 
            # mypose = pyrosetta.pose_from_pdb(f'../out/{protein}/mutated_back_alt.pdb')
            origpose = pyrosetta.pose_from_pdb(f'../out/{s[2]}/{s[1]}')

            stri = mypose.sequence()
            start, end=0,len(stri)
            if s[2] == 'T0974s1':
                start, end = 1, 70
            elif s[2] == 'T1006':
                start, end = 0, 78
            elif s[2] == 'T1072s2':
                start, end = 2,67
            elif s[2] == 'T1046s1':
                start, end = 0, 72
            stri = s[0][start:end]
            
            print('ORIG: ',origpose.sequence())
            print('SEQ : ', mypose.sequence()[start:end])

            ind = origpose.sequence().index(stri)
            print(ind)

            # line up original with my shorter pose
            if ind>0:
                origpose.delete_residue_range_slow(1,ind)
            if end < len(mypose.sequence()):
                mypose.delete_residue_range_slow(end, len(mypose.sequence()))
            if start > 0:
                mypose.delete_residue_range_slow(1,start)
            
            print(f'[log] ######################################################\n origpose sequence: {origpose.sequence()[start: end]}, origpose name: {s[2]}, compare to: {stri[start: end]}, length: {len(stri)}\n\n\n\n')

            print('ORIG: ',origpose.sequence()[0:end-start-1])
            print('SEQ : ', mypose.sequence()[0:end-start-1])
            rmsd = pyrosetta.rosetta.core.scoring.CA_rmsd(mypose, origpose, 0, end-start-1)
            print(f'RMSD: {rmsd}')

            rmsds.append(rmsd)
        w.writerow(rmsds)
            
import csv
import pandas as pd
import os


proteins = [    
                    'T0955' 
                    # 'T1006',
                    # 'T0884', 
                    # 'T0953s1', 
                    # 'T0974s1', 
                    # 'T0980s2', 
                    # 'T1008', 
                    # 'T1046s1', 
                    # 'T1059', 
                    # 'T1072s2'
                ]

# read GD
num_runs = 10
all_min_energies = {}
for p in proteins:        
    min_energies=[]
    for run in range(num_runs):
        filename = '../out/gd_iter70/'+p+'/'+str(run)+'/'+'info.csv'
        if os.path.exists(filename): 
            df_info = pd.read_csv(filename)
            min_energies.append(df_info.iloc[1]['min_energy'])
        else:
            min_energies.append(0)
    min_energies.append(min(min_energies))
    all_min_energies['gd_'+p] = min_energies

# read GD with local MH runs
num_runs = 10
for p in proteins:        
    min_energies=[]
    for run in range(num_runs): 
        filename = '../out/gd_w_localized_65_5/'+p+'/'+str(run)+'/info.csv'
        if os.path.exists(filename):
            df_info = pd.read_csv(filename)
            min_energies.append(df_info.iloc[1]['min_energy'])
        else:
            min_energies.append(0)
    min_energies.append(min(min_energies))
    all_min_energies['gd+localmh_' + p] = min_energies




df = pd.DataFrame(all_min_energies)
df.to_csv('../out/fig2_min_energies_summary.csv')





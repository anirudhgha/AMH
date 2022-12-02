# collect the score histories and minimum energies accross multiple runs to visualize. 

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
#####################################################################################################################################################
# User Input ######################################################################
#####################################################################################################################################################
proteins = [    
                    'T0955',
                    'T0980s2',
                    'T1072s2', 
                    'T0953s1',
                    'T0974s1', 
                    'T1046s1', 
                    'T0884', 
                    'T1006',
                    'T1008' 
                ]
lengths = [
            41, 
            52, 
            69, 
            72, 
            72, 
            74, 
            75, 
            79, 
            80
            ]
folder_names = [    
                    '../out/fig1_mh',
                    '../out/fig1_mh_noisyrestart',
                    '../out/fig1_altmh',
                    '../out/fig1_gd_iter100' ,
                    '../out/fig1_gd_lmh_90_10',
                    '../out/fig1_gd_nr_mhcriterion',
                ]

labels = ['mh', 'mh+nr', 'mh+locmh', 'gd', 'gd+locmh','gd+mhcrit']
colors = ['#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#e6f598', '#abdda4','#66c2a5','#3288bd','#5e4fa2']
# colors_proteins = ['#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#e6f598', '#abdda4','#66c2a5','#3288bd','#5e4fa2']
colors_proteins = ['#feafc7', '#75d6ec', '#010101', '#3b8072', '#ff1d25', '#fe941d', '#D3D303', '#06bd02', '#011a98', '#760188']
skip_methods = []
# plot_proteins = [0]

build_lineplot = False

build_swarmplot = True
boxplot_overlay = False


num_scores_per_history = 10000
num_runs = 10

#####################################################################################################################################################
# Read Data ######################################################################
#####################################################################################################################################################
num_methods = len(folder_names)
samples = {'run': [], 'min_energy':[], 'method':[], 'protein':[], 'length':[]}
sh = [] # sh[method][protein] = (10 x score_history)
for i in range(num_methods):
    if i in skip_methods:
        continue
    count = 0
    sh_per_method = []
    ytemp = np.zeros((1, num_scores_per_history))
    mintemp, protein_per_sample = [], []
    for ip, p in enumerate(proteins):   
        print('processing: ', p)      
        for run in range(num_runs):
            filename = folder_names[i]+'/'+p+'/'+str(run)+'/'+'info.csv'
            if os.path.exists(filename): 
                df_info = pd.read_csv(filename)
                score_history = df_info.loc[:]['score_history'].tolist()[0:num_scores_per_history]
                if np.sum(ytemp) == 0: 
                    ytemp = np.array(score_history)
                else:
                    ytemp = np.vstack((ytemp, np.array(score_history)))

            # store results for this run
            samples['run'].append(count)
            samples['min_energy'].append(np.amin(np.array(score_history)))
            samples['method'].append(labels[i])
            samples['protein'].append(p)
            samples['length'].append(lengths[ip])
            count += 1

        sh_per_method.append(ytemp)
    sh.append(sh_per_method)
sh = np.array(sh)
print(sh.shape)
# print('minimum in T0955 by mhcrit: ', np.min(sh[4]))

#####################################################################################################################################################
# Plot Data 
#####################################################################################################################################################
# mpl.rc('font', family='serif', serif='Times New Roman')
# # Plot 1 : shows energy vs iterations for various MH schemes acting on T0955. 
if build_lineplot:

    fig, ax = plt.subplots(figsize=(10,5))
    # random.seed(9)
    next_method = 0
    x = np.arange(0,num_scores_per_history)
    for i in range(num_methods):
        if i in skip_methods:
                continue
        for p in range(len(proteins)):
            

            # color = ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])][0]
            stds = np.std(sh[next_method][p], axis=0)
            means = np.mean(sh[next_method][p],axis=0)
            minus_std = means-stds
            plus_std = means+stds

            ax.plot(x, means, color=colors[i], alpha=1, label=labels[i], linewidth=1)
            ax.plot(x,minus_std, color=colors[i], alpha=0.15)
            ax.plot(x,plus_std, color=colors[i], alpha=0.15)
            ax.fill_between(x, minus_std, plus_std, color=colors[i], alpha=.1)
        next_method += 1
    ax.grid(b=True, which='major', axis='both', linestyle='-', alpha=0.8)
    # ax.grid(b=True, which='minor', axis='both', linestyle='--', alpha=0.3)

    plottitle= str(proteins)
    ax.set_title('Folding of ' + plottitle)
    ax.set_xlabel('Iterations', fontsize=12)
    ax.set_ylabel('Energy', fontsize=12)
    plt.tick_params(labelsize=12)
    # plt.rcParams.update({'font.size': 14})
    plt.minorticks_on()

    leg = ax.legend()
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    plt.show()

# Plot 2 : Compares the minimum energies of all schemes
if build_swarmplot:
    if len(proteins) == 1:
        fig, ax = plt.subplots(figsize=(5,4))
    else:
        fig, ax = plt.subplots(figsize=(10,5))

    samples = pd.DataFrame(samples)
    if boxplot_overlay:
        ax = sns.boxplot(x='method', y='min_energy', data=samples, showmeans=True, 
                    meanprops= {
                        "marker":"o",
                        "markerfacecolor":"white", 
                        "markeredgecolor":"black"
                        }
                    )
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .1))

    n_colors= len(samples['protein'].unique())
    # palette = sns.color_palette("flare", n_colors=n_colors)

    # append length to protein to make the legend show length
    new_col = []
    custom_palette = {}
    n=0
    for i in range(len(samples['protein'])):
        new_col.append(samples['protein'].iloc[i] + ' - ' + str(samples['length'].iloc[i]))
        if new_col[i] not in custom_palette:
            custom_palette[new_col[i]] = colors_proteins[n]
            n+=1
    samples['legend'] = new_col

    ax = sns.swarmplot(x='method', y='min_energy', hue='legend', palette=custom_palette, data=samples)
    plt.legend(title='Protein', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, frameon=False)

    ax.set_ylabel('Minimum Energy')

    plt.title('Minimization Schemes for ' + str(proteins))
    sns.despine(offset=15)
    # sns.set_theme()
    plt.grid(b=True, which='major', axis='y', linestyle='--', alpha=1, color='k', linewidth=0.5)
    plt.setp(ax.get_xticklabels(), rotation=45)
    # plt.minorticks_on()
    # plt.tick_params(labelsize=12)
    # plt.rcParams.update({'font.size': 14})
    plt.tight_layout()
    plt.show() 

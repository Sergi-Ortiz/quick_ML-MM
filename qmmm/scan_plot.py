#!/Users/sergiortizropero/miniconda3/envs/general/bin/python

# USAGE:
# conda activate general
# python scan_plot.py

import os
import subprocess
bash = '/opt/homebrew/bin/bash'
import numpy as np

# plotting
import matplotlib.pyplot as plt
plt.style.use('default')
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('lines', lw=1, color='b')
rc('legend', loc='best')
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.edgecolor"] = 'black'
plt.rcParams['legend.borderpad'] = 0.25
plt.rcParams['legend.fontsize'] = 11
plt.rcParams.update({'pgf.preamble': r'\usepackage{amsmath}'})

def execute(command_str):
    subprocess.run(command_str, shell=True, executable=bash)

def plot_scan(csv_list, single_csv=True):
    '''
    Plot energy vs rc. 

    Multiple rcs can be plotted. csv files are supplied as [csv1, csv2, etc]
    '''

    # color lists
    colors = ['darkviolet', 'crimson', 'deeppink', 'darkgreen', 'mediumseagreen', 'royalblue', 'darkblue']

    data_list = [np.genfromtxt(csv, delimiter=',') for csv in csv_list]
    
    for i in range(len(csv_list)):

        if single_csv:
            
            # generate time figure
            ncols=2
            fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(6*ncols, 5))

            data = data_list[i]
            _rc = data[1:, 1]
            _energy_hartree = data[1:, 2] / 627.509
            _delta_E = data[1:, 3]
            max_index = np.argmax(_energy_hartree)
            print(f'Highest energy point:\t{max_index+1}')
            print(f'QM/MM energy:        \t{_energy_hartree[max_index]:.6f} Hartree')
            print(r'\D'+f'elta E:            \t{_delta_E[max_index]:.3f} kcal/mol')
            _ch = data[1:, 4]
            _oh = data[1:, 5]
            #print(_delta_E)
            #print(_ch)
            #print(_oh)

            for j in range(ncols):
                if j==0:
                    axes[j].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
                    axes[j].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")
                    axes[j].plot(_rc, _delta_E, color=colors[i], alpha=0.9)
                    
                    axes[j].set(xlabel=r'$\textnormal{Reaction coordinate}\;\Delta r = r_{\textnormal{CH}}-r_{\textnormal{HO}}$', ylabel=r'$\Delta E\;\;(\textnormal{kcal mol}^{-1})$')

                if j==1:
                    axes[j].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
                    axes[j].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")
                     
                    axes[j].plot(_rc, _ch, color='darkgreen', alpha=0.9, label=r'$\textnormal{'+f'C--H'+r' }$')
                    axes[j].plot(_rc, _oh, color='deeppink', alpha=0.9, label=r'$\textnormal{'+f'H--O'+r' }$')
                    
                    axes[j].set(xlabel=r'$\textnormal{Reaction coordinate}\;\Delta r = r_{\textnormal{CH}}-r_{\textnormal{HO}}$', ylabel=r'$\textnormal{Distance (Å)}$')
                    axes[j].legend(loc='best', prop={'size': 10}, title=r'$\textsf{'+f'Distance'+r'}$', title_fontsize=11)
            
        else:
            print('multiple scans not yet suported')
            # TODO plot energy or a SINGLE METRIC (_rc, _delta_E, _ch or _oh) for multiple scans

    plt.tight_layout()
    plt.savefig('./scan_energies.png', dpi=400)




csv_list = ['PES_scan_rc_<FRAME>.csv']
print()
plot_scan(csv_list)
print()
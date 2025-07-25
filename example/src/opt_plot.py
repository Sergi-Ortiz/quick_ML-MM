#!/Users/sergiortizropero/miniconda3/envs/general/bin/python

# USAGE:
# conda activate general
# python opt_plot.py


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

def plot_opt(energies, frames):
    '''
    Plot either just one optimization evolution via energies = [[1., 1.1, 1.05, ...]]
    or multiple runs via energies = [[data_run_1], [data_run_2], ...]

    frames is a list of strings with the frame labels of each optimization run ['frame1', 'frame2', ...]
    '''

    # generate time figure
    fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(6, 5))

    # color lists
    colors = ['darkviolet', 'crimson', 'deeppink', 'darkgreen', 'mediumseagreen', 'royalblue', 'darkblue']

    for i in range(len(energies)):

        axes.tick_params(bottom=True, top=True, left=True, right=False,direction="in")
        axes.tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

        axes.plot(energies[i], color=colors[i], alpha=0.9, label=r'$\textnormal{'+f'{frames[i]}'+r' }$')

    axes.set(xlabel=r'$\textnormal{Optimization step (adim.)}$', ylabel=r'$\textnormal{Energy }\;E\;\;\textnormal{(Hartree)}$')
    axes.legend(loc='best', prop={'size': 10}, title=r'$\textsf{'+f'Frame'+r'}$', title_fontsize=11)

    plt.tight_layout()
    plt.savefig('./opt_energies.png', dpi=400)

execute('grep QM/MM opt.out | awk \'{print $3}\' > opt.txt')
opt_data = [np.loadtxt('opt.txt')]

plot_opt(opt_data, [f'{os.path.split(os.path.split(os.getcwd())[0])[-1]}'])
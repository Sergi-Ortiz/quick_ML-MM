#!/Users/sergiortizropero/miniconda3/envs/mlmm/bin/python
# Sergi Ortiz @ UAB
# 27 July 2025
# ML/MM optimization of a precatalytic structure 

# usage:
# python opt.py -f <FRAME_NUMBER>


#
#   IMPORT SECTION
#
import os
import argparse
import numpy as np

from ase.io import write, read
from ase.io import Trajectory
from ase.visualize import view
from ase.optimize import LBFGS
from ase.constraints import FixAtoms

from mlmm import mlmm_ase

# ------------------------------------------------------------------------

#
#   CUSTOM FUNCTIONS AND CLASSES
#
def freeze_atoms(atoms, freeze_center, freeze_radius, OPT_DIR):

    # transform index to positions
    positions = atoms.positions
    center = positions[freeze_center]     # freeze center is given as in the PDB. Transform to 0-based index

    frozen_idx = [atom.index for atom in atoms if np.linalg.norm(np.array(center) - np.array(positions[atom.index])) > freeze_radius]

    frozen_idx_dir = os.path.join(OPT_DIR, 'frozen_atoms.txt')
    with open(frozen_idx_dir, 'w') as f:
        f.write(f'freeze_atoms: {frozen_idx}')

    # return index list of atoms further than the cutoff
    return frozen_idx

# ------------------------------------------------------------------------


# 
#   CLI ARGUMENTS
# 
parser = argparse.ArgumentParser(description='ML/MM optimization')
parser.add_argument(
    '-f', '--frame', 
    type=str,
    required=True,
    help='Precatalytic frame number.',
)
args = parser.parse_args()


# ------------------------------------------------------------------------


#
#   SIMULATION PARAMETERS
#

# frame src
frame = args.frame
SRC_PATH = os.path.join('../', f'{frame}', 'src')            # frame/src
OPT_PATH = os.path.join('../', f'{frame}', 'opt')            # frame/opt

pdb_path = os.path.join(SRC_PATH, f'{frame}_solv_cropped.pdb')           # .pdb from where the prmtop and inpcrd where generated (frame src)
prmtop_path = os.path.join(SRC_PATH, f'{frame}_solv_cropped.prmtop')     # .prmtop
inpcrd_path = os.path.join(SRC_PATH, f'{frame}_solv_cropped.inpcrd')     # .inpcrd (or rst7 restart file???)

# opt files
opt_log = os.path.join(OPT_PATH, f'opt_{frame}.log')
opt_pdb = os.path.join(OPT_PATH, f'opt_{frame}.pdb')        # resulting optimized structure path


# source ML/MM region
ml_region_path = os.path.join(SRC_PATH, 'ml_region.pdb')                 # .pdb of the ML region


# THEORY OPTIONS
ml_charge = 1           # charge of ML region including link atoms
ml_spin = 6             # multiplicity of ML region
mm_threads = 12          # CPU threads to use for the MM region

# OPT OPTIONS
opt_steps = 2000
opt_thereshold = 0.001   # eV / Å
freeze_radius = 15.0     # radius from the H from which to freeze MM atoms (Å). SHOULD BE 15

# C_idx = 10581         # C10
# C_idx = 10574         # C13
freeze_center = 10581-1   # index of the C_abs atom to serve as center (same as in QM/MM). 1-BASED INDEXING!!!

# ------------------------------------------------------------------------

# starting opt structure (the same as the cropped one)
initial_pdb = pdb_path
#initial_pdb = os.path.join(SRC_PATH, f'scan_{frame}_14620.pdb')
#initial_pdb = os.path.join(SRC_PATH, 'mod_TS.pdb')

# generate scan path
if not os.path.exists(OPT_PATH):
    os.makedirs(OPT_PATH)

# ML/MM calculator setup
mlmm_kwargs = dict(
    real_pdb     = pdb_path,
    real_parm7   = prmtop_path,
    real_rst7    = inpcrd_path,
    model_pdb    = ml_region_path,
    model_charge = ml_charge,       
    model_mult   = ml_spin,         
    backend      = 'uma',           # "uma" or "aimnet2"
    uma_model    = 'uma-s-1p1',
    ml_device    = 'auto',          # "auto" | "cuda" | "cpu"
    ml_cuda_idx  = 0,
    mm_device    = 'cuda',
    mm_cuda_idx  = 0,
    mm_threads   = mm_threads,      # max number of threads available
)

atoms = read(initial_pdb)
atoms.calc = mlmm_ase(**mlmm_kwargs)

# freeze far away atoms
frozen_idx = freeze_atoms(atoms, freeze_center, freeze_radius, OPT_PATH)
freeze_c = FixAtoms(indices=frozen_idx)
atoms.set_constraint(freeze_c)

# optimization
opt = LBFGS(atoms, logfile=opt_log)
opt.run(fmax=opt_thereshold, steps=opt_steps)

# save to file
write(opt_pdb, atoms)
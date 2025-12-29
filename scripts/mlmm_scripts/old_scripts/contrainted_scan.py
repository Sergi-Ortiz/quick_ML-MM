#!/Users/sergiortizropero/miniconda3/envs/mlmm/bin/python
# Sergi Ortiz @ UAB
# 27 July 2025
# Perform constrained ML/MM scan

# usage:
# python constrained_scan.py -f <FRAME_NUMBER>


#
#   IMPORT SECTION
#
import os
import argparse
import numpy as np

from ase.io import write, read
from ase.optimize import LBFGS
from ase.constraints import FixAtoms

from mlmm import mlmm_ase

# ------------------------------------------------------------------

#
#   CUSTOM FUNCTIONS AND CLASSES
#
class HarmonicConstraint2D:
    """Explore a reaction coordinate based on a 1D distance with harmonic constraints"""
        
    def __init__(self, atoms, a_H, a_C, a_O, k, step):

        '''
        # model H tranfer between C and O (C-H  O --> C  H-O)
        a_H position of atom 1 (H, 'recipient' of the contraint)
        a_C position of atom 2 (C bonded atom, bonded to H at the start) 
        a_O position of atom 3 (O target atom, bonded to H at the end)
        k harmonic constant eV / Å^2
        step diference in reaction coordinate 
        '''

        # indices of the atoms
        self.a_H = a_H
        self.a_C = a_C
        self.a_O = a_O

        # initial positions when generating the constraint
        self.positions_0 = atoms.positions
        self.pos_H_0 = self.positions_0[self.a_H]
        self.pos_C_0 = self.positions_0[self.a_C]
        self.pos_O_0 = self.positions_0[self.a_O]

        # force constant
        self.k = k
        self.step = step

        # reaction coordinate
        self.r_CH_0 = np.linalg.norm(np.array(self.pos_H_0) - np.array(self.pos_C_0))
        self.r_OH_0 = np.linalg.norm(np.array(self.pos_H_0) - np.array(self.pos_O_0))
        self.r_0 = self.r_CH_0 - self.r_OH_0


    def adjust_positions(self, atoms, newpositions):
        pass

    def adjust_momenta(self, atoms, momenta):
        pass


    def adjust_potential_energy(self, atoms):
        positions = atoms.positions
        pos_H = positions[self.a_H]
        pos_C = positions[self.a_C]
        pos_O = positions[self.a_O]

        r_CH = np.linalg.norm(np.array(pos_H) - np.array(pos_C))
        r_OH = np.linalg.norm(np.array(pos_H) - np.array(pos_O))

        E = self.k/2. * ((r_CH - r_OH) - (self.r_CH_0 - self.r_OH_0) - self.step)**2
        return E


    def adjust_forces(self, atoms, forces):
        positions = atoms.positions
        pos_H = positions[self.a_H]
        pos_C = positions[self.a_C]
        pos_O = positions[self.a_O]

        r_CH = np.linalg.norm(np.array(pos_H) - np.array(pos_C))
        r_OH = np.linalg.norm(np.array(pos_H) - np.array(pos_O))

        # return the gradient of the energy
        F_prefactor = - self.k * ((r_CH - self.r_CH_0) - (r_OH - self.r_OH_0) - self.step) 
        F_vector = ((pos_H - pos_C)/(r_CH) - (pos_H - pos_O)/(r_OH))
        F = F_prefactor * F_vector

        # three components
        forces[self.a_H] += F 


def perform_step(atoms, l, frame, SCAN_DIR, H_idx, C_idx, O_idx, k, rc_step, opt_thereshold=0.05, opt_steps=10000, frozen_idx=[]):
    '''
    # SCAN SINGLE STEP (CONSTRAINED OPTIMIZATION)
    # harmonic constraint to advance in the reaction coordinate
    # perform a local optimization (atoms beyond a certain thereshold are frozen)
    # freeze atoms for the optimization
    '''

    # suppose input is 0-based indexing!!!!


    # transform 1-based indexing to 0-based for the constraint
    harmonic_c = HarmonicConstraint2D(atoms=atoms, a_H=H_idx, a_C=C_idx, a_O=O_idx, k=k, step=rc_step)

    if len(frozen_idx) != 0:
        freeze_c = FixAtoms(indices=frozen_idx)

    # freeze Fe atom
    Fe_idx = [atom.index for atom in atoms if atom.symbol == 'Fe']
    #print(Fe_idx)
    Fe_c = FixAtoms(indices=Fe_idx)

    # freeze OH1 O atom (ad hoc)
    OH_c = FixAtoms(indices=[O_idx])

    # combine constraints
    atoms.set_constraint()
    atoms.set_constraint([harmonic_c, freeze_c, OH_c, Fe_c])

    log_k_path = os.path.join(SCAN_DIR, f'scan_{frame}_{l}.log')
    opt = LBFGS(atoms, logfile=log_k_path)
    opt.run(fmax=opt_thereshold, steps=opt_steps)

    energy = atoms.get_potential_energy()

    return atoms, energy
    


def compute_rc(atoms, H_idx, C_idx, O_idx):
    '''
    computes the reaction coordinate of a given geometry
    '''

    # transform 1-based indexing to 0-based
    positions = atoms.positions
    pos_H = positions[H_idx]
    pos_C = positions[C_idx]
    pos_O = positions[O_idx]

    r_CH = np.linalg.norm(np.array(pos_H) - np.array(pos_C))
    r_OH = np.linalg.norm(np.array(pos_H) - np.array(pos_O))
    _rc = r_CH - r_OH

    return _rc, r_CH, r_OH


def freeze_atoms(atoms, freeze_center, freeze_radius):

    # transform 1-based indexing to 0-based
    positions = atoms.positions
    center = positions[freeze_center]     # freeze center is given as in the PDB. Transform to 0-based index

    # return index list of atoms further than the cutoff
    return [atom.index for atom in atoms if np.linalg.norm(np.array(center) - np.array(positions[atom.index])) > freeze_radius]

# ------------------------------------------------------------------


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
SCAN_PATH = os.path.join('../', f'{frame}', 'scan')          # frame/scan
OPT_PATH = os.path.join('../', f'{frame}', 'opt')

#pdb_path = os.path.join(SRC_PATH, f'scan_{frame}_54.pdb')               # restart scan from a given pdb
pdb_path = os.path.join(OPT_PATH, f'opt_{frame}.pdb')                    # .pdb from where the prmtop and inpcrd where generated (frame src)
prmtop_path = os.path.join(SRC_PATH, f'{frame}_solv_cropped.prmtop')     # .prmtop
inpcrd_path = os.path.join(SRC_PATH, f'{frame}_solv_cropped.inpcrd')     # .inpcrd (or rst7 restart file???)
csv_out = os.path.join(SCAN_PATH, f'scan_{frame}.csv')                   # output .csv file

# source ML/MM region
ml_region_path = os.path.join(SRC_PATH, 'ml_region.pdb')                 # .pdb of the ML region

# THEORY OPTIONS
ml_charge = 1                       # charge of ML region including link atoms
ml_spin = 6                         # multiplicity of ML region
mm_threads = 12                     # CPU threads to use for the MM region

# other opt options
opt_steps = 250                     # max opt steps
opt_thereshold = 0.08               # eV / Å
freeze_radius = 15.0                # radius from the H from which to freeze MM atoms (Å). SHOULD BE 15

# constraints options (1 BASED INDEXING!) (from the PDB)
# C10 abstraction
H_idx = 10583-1         # index of H atom
C_idx = 10581-1         # index of start bond atom X-H     Y
O_idx = 10552-1         # index of final bond atom X     H-Y 

# C13 abstraction
#H_idx = 10576         # index of H atom
#C_idx = 10574         # index of start bond atom X-H     Y
#O_idx = 10552         # index of final bond atom X     H-Y 

k = 291.54            # eV/Å^2
#k = 300.0

# scan options
rc_step = 0.01               # Å
rc_thereshold = 3.0         # stop scan when rc thereshold is reached
visualize = False

# starting opt structure 
opt_pdb = pdb_path

# generate scan path
if not os.path.exists(SCAN_PATH):
    os.mkdir(SCAN_PATH)


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
    mm_device    = 'cpu',
    mm_cuda_idx  = 0,
    mm_threads   = mm_threads,      # max number of threads available
)
CALC_ASE = mlmm_ase(**mlmm_kwargs)


# ------------------------------------------------------------------

#
#   ACTUAL SIMULATION
#
atoms = read(opt_pdb)
atoms.calc = CALC_ASE

with open(csv_out, 'w') as f:
    f.write('j,rc [Å],energy [eV],r_CH,r_OH\n')

# obtain initial rc
_rc, r_CH, r_OH = compute_rc(atoms, H_idx, C_idx, O_idx)

# iterate all reaction coordinates
l = 1
while _rc < rc_thereshold:

    

    # get frozen indices
    frozen_idx = freeze_atoms(atoms, H_idx, freeze_radius)

    # perform the optimization
    atoms, energy = perform_step(
        atoms, 
        l, 
        frame, 
        SCAN_PATH,
        H_idx, 
        C_idx, 
        O_idx, 
        k, 
        rc_step, 
        opt_thereshold=opt_thereshold, 
        opt_steps=opt_steps, 
        frozen_idx=frozen_idx
    )

    # save to file
    new_structure = os.path.join(SCAN_PATH, f'scan_{frame}_{l}.pdb')
    write(new_structure, atoms)

    with open(csv_out, 'a') as f:
        f.write(f'{l},{_rc:.6f},{energy:.6f},{r_CH},{r_OH}\n')
    
    #atoms = read(new_structure)
    #atoms.calc = CALC_ASE

    # obtain rc of the new structure
    _rc, r_CH, r_OH = compute_rc(atoms, H_idx, C_idx, O_idx)

    # save results to a .csv (same as in QM/MM)
    
    print(f'Finished step {l}\n')
    l += 1

else:
    print(f'Constrained scan ended normally at {_rc}')



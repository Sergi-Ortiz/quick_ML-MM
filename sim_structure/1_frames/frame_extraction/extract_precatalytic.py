#!/Users/sergiortizropero/miniconda3/envs/md_traj/bin/python
# Sergi Ortiz @ UAB
# 25 July 2025
# Extract precatalytic structures from MD trajectory

# example usecase:
# python extract_precatalytic.py -m 5-HETE -H 10 -d 4.0 -l 20.0 -n 1 


# 
#   IMPORT SECTION
#
import os
import argparse
import numpy as np

# md analysis
from MDAnalysis import Universe
from EMDA import EMDA, __version__
from EMDA.frame_extractor import trajectory_frame_extractor
print(f'EMDA version:\t{__version__}')


# 
#   CLI ARGUMENTS
# 
parser = argparse.ArgumentParser(description='Extract precatalytic frames from a MD trajectory')

# positional arguments
parser.add_argument(
    '-H', '--hydrogen', 
    type=int,
    required=True,
    help='Hydrogen selection for which to perform the QM/MM scan. Options: 10, 13.',
)
parser.add_argument(
    '-d', '--distance',
    default=4.0,
    type=float,
    required=False,
    help='Distance criterium (H--O cofactor) to use for the precatalytic definition. Default 4.0',
)
parser.add_argument(
    '-l', '--dihedral',
    default=20.0,
    type=float,
    required=False,
    help='Dihedral criterium (around H atom) to use for the precatalytic definition, plus minus in deg. Default 20.0',
)
parser.add_argument(
    '-m', '--molecule', 
    type=str,
    required=True,
    help='Ligand used in the MD trajectory. Options: \'AA\', \'5-HETE\', \'5-HpETE\'.',
)
parser.add_argument(
    '-n', '--number', 
    type=int,
    default=1,
    required=False,
    help='Number of precatalytic structures to extract from the MD simulation. The selection is random. Default: 1.',
)
parser.add_argument(
    '-v', '--verbose', 
    type=bool,
    default=False,
    required=False,
    help='Whether to verbose. Default: False',
)

# obtain the arguments
args = parser.parse_args()
selection = args.hydrogen
dist_criterium = args.distance
dih_criterium = args.dihedral
lig_ID = args.molecule
num_precat = args.number
verbose = args.verbose


# constants
md_base_dir = os.path.abspath('./')
parent_dir = 'md_run'   
protein_ID = 'hLOX15'       


#
#   EXTRACT THE MD DATA
#
lig_ID_list = ['AA', '5-HpETE', '5-HETE']
if lig_ID not in lig_ID_list:
    raise ValueError(f'Introduced lig_ID {lig_ID} is not among the supported ligands: {lig_ID_list}')

# topology
topologies_1 = [i for i in os.listdir(os.path.join(md_base_dir, parent_dir)) if i.endswith('.prmtop')]
top = os.path.join(md_base_dir, parent_dir, topologies_1[0])
if len(topologies_1) > 1:
    raise ValueError('bruh. More than one topology files? Check md source directory.')
print(f'loading prmtop\t{top}')

# production trajectory
prod_dir_1 = os.path.join(md_base_dir, parent_dir, 'prod/out')                    # locate the production traj
prod_files = [int(x.strip('.nc').strip('prod_')) for x in os.listdir(prod_dir_1) if '.nc' in x]
traj_1 = [os.path.join(prod_dir_1, f'prod_{i}.nc')  for i in prod_files]      # assume constant length 250 ns, 25 
print(np.array(traj_1))

# load base MD (r=1)
emda = EMDA()
emda.load_variant(
    parameters=top,
    trajectory=traj_1,
    variant_name=protein_ID
)
u = Universe(top, traj_1)
total_time_simulated = 10 * len(traj_1)
print(f'total time simulated: {total_time_simulated} ns')


#
#   SELECTIONS
#
print(f'performing atom selections for {lig_ID}')

# AA
if lig_ID == lig_ID_list[0]:
    if selection == 10:
        emda.select('H10A', 'name H10A and resid 665')      # pro-S H10
        #emda.select('H10B', 'name H10B and resid 665')     # pro-R H10

        emda.select('p10_A', 'name C12 and resid 665')      # A i B formen el cis double bond
        emda.select('p10_B', 'name C11 and resid 665')
        emda.select('p10_C', 'name C10 and resid 665')      # C conté H10A i H10B
        emda.select('p10_D', 'name C9 and resid 665')       # D és el carboni més proxim al COO

    elif selection == 13:
        emda.select('H13A', 'name H13A and resid 665')      # pro-S H13
        #emda.select('H13B', 'name H13B and resid 665')     # pro-S H13

        emda.select('p13_A', 'name C11 and resid 665')      # A i B formen el cis double bond
        emda.select('p13_B', 'name C12 and resid 665')
        emda.select('p13_C', 'name C13 and resid 665')      # C conté H13A i H13B
        emda.select('p13_D', 'name C14 and resid 665')      # D és el carboni més llunyà al COO

# 5-HpETE
elif lig_ID == lig_ID_list[1]:
    if selection == 10:
        emda.select('H10A', 'name H19 and resid 665')      # pro-S H10
        #emda.select('H10B', 'name H18 and resid 665')     # pro-R H10

        emda.select('p10_A', 'name C9 and resid 665')       # A i B formen el cis double bond
        emda.select('p10_B', 'name C10 and resid 665')    
        emda.select('p10_C', 'name C11 and resid 665')      # C conté H10A i H10B
        emda.select('p10_D', 'name C12 and resid 665')      # D és el carboni més proxim al COO

    elif selection == 13:
        emda.select('H13A', 'name H15 and resid 665')      # pro-S H13
        #emda.select('H13B', 'name H14 and resid 665')     # pro-S H13

        emda.select('p13_A', 'name C10 and resid 665')      # A i B formen el cis double bond
        emda.select('p13_B', 'name C9 and resid 665')
        emda.select('p13_C', 'name C8 and resid 665')       # C conté H13A i H13B
        emda.select('p13_D', 'name C7 and resid 665')       # D és el carboni més llunyà al COO

# 5-HETE
elif lig_ID == lig_ID_list[2]:


    if selection == 10:
        emda.select('H10A', 'name H19 and resid 665')      # pro-S H10
        #emda.select('H10B', 'name H18 and resid 665')     # pro-R H10

        emda.select('p10_A', 'name C9 and resid 665')       # A i B formen el cis double bond
        emda.select('p10_B', 'name C10 and resid 665')    
        emda.select('p10_C', 'name C11 and resid 665')      # C conté H10A i H10B
        emda.select('p10_D', 'name C12 and resid 665')      # D és el carboni més proxim al COO

    elif selection == 13:
        emda.select('H13A', 'name H15 and resid 665')      # pro-S H13
        #emda.select('H13B', 'name H14 and resid 665')     # pro-S H13

        emda.select('p13_A', 'name C10 and resid 665')      # A i B formen el cis double bond
        emda.select('p13_B', 'name C9 and resid 665')
        emda.select('p13_C', 'name C8 and resid 665')       # C conté H13A i H13B
        emda.select('p13_D', 'name C7 and resid 665')       # D és el carboni més llunyà al COO

else:
    raise ValueError('Something went wrong selecting the atoms to be analyzed')

emda.select('OH1', 10552, sel_type='at_num')        # O atom from the hydroxyl cofactor
print('\n\nselected selections:\n')
print(emda.selections)

# define H10 distance and planarity
if selection == 10:

    emda.add_distance('distance H10A -- O OH1', 'OH1', 'H10A')
    #emda.add_distance('distance H10B -- O OH1', 'OH1', 'H10B')

    emda.add_dihedral('dihedral C12--C11--C10--C9 (planarity around C10)', 'p10_A', 'p10_B', 'p10_C', 'p10_D', domain=360)  

# define H13 distance and planarity
elif selection == 13:

    emda.add_distance('distance H13A -- O OH1', 'OH1', 'H13A')
    #emda.add_distance('distance H13B -- O OH1', 'OH1', 'H13B')

    emda.add_dihedral('dihedral C11--C12--C13--C14 (planarity around C13)', 'p13_A', 'p13_B', 'p13_C', 'p13_D', domain=360)

print('\n\nselected measurements:\n')
print(emda.measures)


#
#   TRAJECTORY PREPROCESSING
#
precat_pkl = os.path.join(md_base_dir, f'{protein_ID}_{lig_ID}_H{selection}_precat_analysis.pkl')
try:
    # if previously computed, load the pkl file  
    print(f'loading analysis from {precat_pkl}')
    emda.load(precat_pkl)
except:
    # analyze the trajectory from scratch 
    print(f'analyzing precatalytic structures in the trajectory {parent_dir}')
    emda.run()
    emda.save(precat_pkl[:-4])    


#
#   PERFORM PRECATALYTIC ANALYSIS
#
if selection == 10:
    print('analyzing C10 dist and dihedral')
    emda.analyse_value('dihedral_criterion', 'dihedral C12--C11--C10--C9 (planarity around C10)', mode='thres', val1=180-dih_criterium, val2=180+dih_criterium)
    emda.analyse_value('distance_criterion', 'distance H10A -- O OH1', mode='thres', val1=dist_criterium)
    
elif selection == 13:
    print('analyzing C13 dist and dihedral')
    emda.analyse_value('dihedral_criterion', 'dihedral C11--C12--C13--C14 (planarity around C13)', mode='thres', val1=180-dih_criterium, val2=180+dih_criterium)
    emda.analyse_value('distance_criterion', 'distance H13A -- O OH1', mode='thres', val1=dist_criterium)

else:
    raise ValueError('ep! Analysis step exploded!')

# perform the analysis
criterion_1 = list(emda.analyses['dihedral_criterion'].result.values())[0]['R1']
criterion_1_idx = np.where(criterion_1)[0]
criterion_2 = list(emda.analyses['distance_criterion'].result.values())[0]['R1']
criterion_2_idx = np.where(criterion_2)[0]

# AND-type restriction
criteria = list(np.array([criterion_1, criterion_2]).min(axis=0))
criteria_idx = np.where(criteria)[0]

#print(f'Dihedral criterion\n{criterion_1}\n')
#print(f'Distance criterion\n{criterion_2}\n')
#print(f'Both ({criteria.count(True)} frames)\n{criteria}\n')

# final results
if verbose:
    print(f'Frames that satisfy the dihedral criterion:\t{(len(criterion_1_idx)*100/(1000 * len(traj_1))):.2f}\t%')
    print(f'Frames that satisfy the distance criterion:\t{(len(criterion_2_idx)*100/(1000 * len(traj_1))):.2f}\t%')
    print(f'Frames that satisfy both criteria:\t\t{(len(criteria_idx)*100/(1000 * len(traj_1))):.2f}\t%')


#
#   FRAME EXTRACTION
#
# draw a random number of precatalytic structures to extract
rng = np.random.default_rng(seed=123)
try:
    random_precat_idx = criteria_idx[rng.choice(len(criteria_idx), size=num_precat, replace=False, )].tolist()
except:
    # number of frames is less than N to be extracted 
    print(f'Unable to extract {num_precat} from {len(criteria_idx)} total precatalytic structures. Extracting {len(criteria_idx)}.')
    random_precat_idx = criteria_idx[rng.choice(len(criteria_idx), size=len(criteria_idx), replace=False)].tolist()
    print(random_precat_idx)

# extract the corresponding selected frames
frame_path = f'./frames_H{selection}'
trajectory_frame_extractor(u, random_precat_idx, folder=frame_path, outformat=".pdb", solvent=True)

print('Frame extraction finished normally.')
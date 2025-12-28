#!/Users/sergiortizropero/miniconda3/envs/md_traj/bin/python
# Sergi Ortiz @ UAB
# 25 July 2025
# Prepare all frames for QM/MM calculations.
# Frames are automatically copied to QM-MM directory for further setup with 
# setup_qmmm.py, which actually sets up the SLURM jobs. 

# example usecase:
# python prepare_frames.py -H 10 -m 5-HETE

# 
#   IMPORT SECTION
#
import os
import shutil
import subprocess
import argparse


#
#   FUNCTIONS
#
bash = '/opt/homebrew/bin/bash'     # macOS
#bash = '/bin/bash'                 # linux

def execute(command_str):
    subprocess.run(command_str, shell=True, executable=bash)


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
    '-m', '--molecule', 
    type=str,
    required=True,
    help='Ligand used in the MD trajectory. Options: \'AA\', \'5-HETE\', \'5-HpETE\'.',
)
args = parser.parse_args()
H = args.hydrogen
lig_ID = args.molecule


#
#   QM-MM INFORMATION
#

# set up correct atoms for the setup depending on ligand and H number
if lig_ID == '5-HETE':
    if H == 10:
        gen_params_atom = 10581
        H_atom = 10583
        C_atom = 10581
        O_OH1_atom = 10552
    elif H == 13:
        gen_params_atom = 10574
        H_atom = 10576
        C_atom = 10574
        O_OH1_atom = 10552
    else:
        raise ValueError(f'H{H} is not H10 or H13!')
    
elif lig_ID == '5-HpETE':
    if H == 10:
        gen_params_atom = 10581
        H_atom = 10583
        C_atom = 10581
        O_OH1_atom = 10552
    elif H == 13:
        gen_params_atom = 10574
        H_atom = 10576
        C_atom = 10574
        O_OH1_atom = 10552
    else:
        raise ValueError(f'H{H} is not H10 or H13!')
else:
    raise ValueError(f'failed to identify ligand name {lig_ID} as 5-HETE or 5-HpETE.')


#
#   PATH SETUP
#
FRAME_DIR = os.path.relpath(f'../frame_extraction/frames_H{H}')
frame_parametrization = '../frame_parametrization'
frame_save_qmmm = '../../QM-MM'
frame_save_mlmm = '../../ML-MM'

# QM-MM specific files
qm_list_name = f'qm_list_{lig_ID[2:]}'
tleap_name = f'tleap_{lig_ID[2:]}'



#
#   GENERATE PARAMETERS FOR EACH PRECATALYTIC STRUCTURE
#

selected_list = [int(x.strip('.pdb')) for x in os.listdir(FRAME_DIR) if ('.pdb' in x) and ('_' not in x)]

# create frame_preprocessing and copy src code
if not os.path.exists(frame_parametrization):
    os.mkdir(frame_parametrization)

# copy necessary preprocessing scripts and tleap folder with params. 
execute(f'cp -r {tleap_name} {os.path.join(frame_parametrization, 'tleap')}')
execute(f'cp amber_to_chemshell_modeller.py {frame_parametrization}')
execute(f'cp generate_params.sh {frame_parametrization}')
execute(f'cd {frame_parametrization}; sed "s/<CATOM>/{C_atom}/g" generate_params.sh > generate_params_mod.sh; rm generate_params.sh; mv generate_params_mod.sh generate_params.sh')

# copy frames to their respective folders
for frame in selected_list:
    execute(f'cp {os.path.join(FRAME_DIR, str(frame)+'.pdb')} {frame_parametrization}')

# process the frames with tleap via `generate_params.sh` and `amber_to_chemshell_modeller.py`
execute(f'cd {frame_parametrization}; bash generate_params.sh')


#
#   QM-MM tidy setup
#
# save .prmtop .inpcrd .pdb and .act_list in a folder for each preprocessed frame
for frame in selected_list:
    if not os.path.exists(os.path.join(frame_save_qmmm, str(frame))):
        os.makedirs(os.path.join(frame_save_qmmm, str(frame), 'src'))
        
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv_cropped_mod.prmtop ../../../../QM-MM/{frame}/src')
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv_cropped.inpcrd ../../../../QM-MM/{frame}/src')
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv_cropped.pdb ../../../../QM-MM/{frame}/src')
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv.act_list ../../../../QM-MM/{frame}/src')

# copy qmmm_scripts to QM-MM folder
shutil.copytree('../../scripts/qmmm_scripts', '../../QM-MM/qmmm_scripts')

#
#   ML-MM tidy setup
#
# save original trimmed .prmtop .inpcrd .pdb to a folder 
for frame in selected_list:
    if not os.path.exists(os.path.join(frame_save_mlmm, str(frame))):
        os.makedirs(os.path.join(frame_save_mlmm, str(frame)))
        os.mkdir(os.path.join(frame_save_mlmm, str(frame), 'src'))
        
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv_cropped.prmtop ../../../../ML-MM/{frame}/src')
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv_cropped.inpcrd ../../../../ML-MM/{frame}/src')
    execute(f'cd {os.path.join(frame_parametrization, str(frame), 'tleap')}; cp {frame}_solv_cropped.pdb ../../../../ML-MM/{frame}/src')

# copy mlmm_scripts to ML-MM folder
shutil.copytree('../../scripts/mlmm_scripts', '../../ML-MM/mlmm_scripts')



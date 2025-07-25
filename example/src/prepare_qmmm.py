#!/Users/sergiortizropero/miniconda3/envs/md_traj/bin/python
# Sergi Ortiz @ UAB
# 25 July 2025
# Extract precatalytic structures from MD trajectory

# example usecase:
# python prepare_qmmm.py -m 5-HETE -H 10 -d 4.0 -l 20.0 -n 1 



# ! ! !    IMPORTANT    ! ! !
# THE PATH TREE FOR THIS SCRIPT IS VERY SPECIAL. REFER TO AN EXAMPLE @ GITHUB
# required directories
#   - 1_frames/frame_extraction/
#       - md_run, frames.ipynb, frames
#   - src
#       - qm/mm setup and source code




# 
#   IMPORT SECTION
#
import os
import subprocess
import argparse


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

# 
#   QM/MM SETUP OPTIONS
#
args = parser.parse_args()
H = args.hydrogen
dist_criterium = args.distance
dih_criterium = args.dihedral
lig_ID = args.molecule
num_precat = args.number



# macOS
bash = '/opt/homebrew/bin/bash'
# linux
#bash = '/bin/bash'

def execute(command_str):
    subprocess.run(command_str, shell=True, executable=bash)


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

# FRAMES
frame_extraction = '../1_frames/frame_extraction'
frame_path = os.path.relpath(f'../1_frames/frame_extraction/frames_H{H}')
frame_preprocessing = '../1_frames/frame_preprocessing'
frame_source = '../1_frames/frame_source'

# CALCULATIONS
calculations = '../2_calculations'


# ligand specific documents
qm_list_name = f'qm_list_{lig_ID[2:]}'
tleap_name = f'tleap_{lig_ID[2:]}'


#
#   PREPROCESS THE MD TRAJECTORY AND EXTRACT PRECATALYTIC STRUCTURES
#
# copy extract_precatalytic.py to frame_extraction and execture the script 
execute(f'cp extract_precatalytic.py {frame_extraction}')
execute(f'cd {frame_extraction}; python extract_precatalytic.py -m {lig_ID} -H {H} -d {dist_criterium} -l {dih_criterium} -n {num_precat} 2>&1 | tee extract_H{H}.txt')


#
#   GENERATE PARAMETERS FOR EACH PRECATALYTIC STRUCTURE
#

# select frames_int random frames
selected_list = [int(x.strip('.pdb')) for x in os.listdir(frame_path) if ('.pdb' in x) and ('_' not in x)]

# create frame_preprocessing and copy src code
if not os.path.exists(frame_preprocessing):
    os.mkdir(frame_preprocessing)

# copy necessary preprocessing scripts and tleap folder with params. 
execute(f'cp -r {tleap_name} {os.path.join(frame_preprocessing, 'tleap')}')
execute(f'cp amber_to_chemshell_modeller.py {frame_preprocessing}')
execute(f'cp generate_params.sh {frame_preprocessing}')
execute(f'cd {frame_preprocessing}; sed "s/<CATOM>/{C_atom}/g" generate_params.sh > generate_params_mod.sh; rm generate_params.sh; mv generate_params_mod.sh generate_params.sh')

# copy frames to their respective folders
for frame in selected_list:
    execute(f'cp {os.path.join(frame_path, str(frame)+'.pdb')} {frame_preprocessing}')

# process the frames with tleap via `generate_params.sh` and `amber_to_chemshell_modeller.py`
execute(f'cd {frame_preprocessing}; bash generate_params.sh')


# add relevant parameters as src for each frame 
if not os.path.exists(frame_source):
    os.mkdir(frame_source)

for frame in selected_list:
    if not os.path.exists(os.path.join(frame_source, str(frame))):
        os.mkdir(os.path.join(frame_source, str(frame)))
        
    execute(f'cd {os.path.join(frame_preprocessing, str(frame), 'tleap')}; cp {frame}_solv_cropped_mod.prmtop ../../../frame_source/{frame}')
    execute(f'cd {os.path.join(frame_preprocessing, str(frame), 'tleap')}; cp {frame}_solv_cropped.inpcrd ../../../frame_source/{frame}')
    execute(f'cd {os.path.join(frame_preprocessing, str(frame), 'tleap')}; cp {frame}_solv_cropped.pdb ../../../frame_source/{frame}')
    execute(f'cd {os.path.join(frame_preprocessing, str(frame), 'tleap')}; cp {frame}_solv.act_list ../../../frame_source/{frame}')


# prepare the QM/MM setup for each frame
for frame in selected_list:
    if not os.path.exists(os.path.join(calculations, str(frame))):
        
        # 1_qm_region
        os.makedirs(os.path.join(calculations, str(frame)))
        os.mkdir(os.path.join(calculations, str(frame), '1_qm_region'))
        execute(f'cp {qm_list_name} {os.path.join(calculations, str(frame), '1_qm_region', 'qm_list')}')
        execute(f'cp {os.path.join(frame_source, str(frame), f'{frame}_solv_cropped.pdb')} {os.path.join(calculations, str(frame), '1_qm_region')}')
        execute(f'cp atoms_list_to_pdb.py {os.path.join(calculations, str(frame), '1_qm_region')}')
        print('converting qm_list to qm_list.pdb')
        execute(f'cd {os.path.join(calculations, str(frame), '1_qm_region')}; python atoms_list_to_pdb.py \'qm_list\' {frame}_solv_cropped.pdb qm_list')

        # open chimera with the #0 qm_list.pdb and #1 the protein
        # run the script
        execute(f'cp qm_region_view.py {os.path.join(calculations, str(frame), '1_qm_region')}')
        execute(f'cd {os.path.join(calculations, str(frame), '1_qm_region')}; sed "s:<PATH>:{os.path.abspath(os.path.join(calculations, str(frame), '1_qm_region', 'view.py'))}:g" qm_region_view.py > qm_region_view_mod.py; rm qm_region_view.py; mv qm_region_view_mod.py qm_region_view.py')

        # 2_opt
        os.mkdir(os.path.join(calculations, str(frame), '2_opt'))
        execute(f'cp {os.path.join(frame_source, str(frame), '*')} {os.path.join(calculations, str(frame), '2_opt')}')
        execute(f'cp {qm_list_name} {os.path.join(calculations, str(frame), '2_opt', 'qm_list')}')
        execute(f'cp opt_job.sh {os.path.join(calculations, str(frame), '2_opt')}')
        execute(f'cp opt_plot.py {os.path.join(calculations, str(frame), '2_opt')}')
        execute(f'cp opt.chm {os.path.join(calculations, str(frame), '2_opt')}')
        execute(f'cd {os.path.join(calculations, str(frame), '2_opt')}; sed "s/<FRAME>/{frame}/g" opt_job.sh > opt_job_mod.sh; rm opt_job.sh; mv opt_job_mod.sh opt_job.sh')
        execute(f'cd {os.path.join(calculations, str(frame), '2_opt')}; sed "s/<FRAME>/{frame}/g" opt.chm > opt_mod.chm; rm opt.chm; mv opt_mod.chm opt.chm')

        # 3_scan
        os.mkdir(os.path.join(calculations, str(frame), '3_scan'))
        execute(f'cp {os.path.join(os.path.join(calculations, str(frame), '2_opt'), f'{frame}_solv_cropped_mod.prmtop')} {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cp {os.path.join(os.path.join(calculations, str(frame), '2_opt'), f'{frame}_solv_cropped.inpcrd')} {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cp {os.path.join(os.path.join(calculations, str(frame), '2_opt'), f'{frame}_solv.act_list')} {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cp {os.path.join(os.path.join(calculations, str(frame), '2_opt'), f'qm_list')} {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cp scan_job.sh {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cp scan_plot.py {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cp scan.chm {os.path.join(calculations, str(frame), '3_scan')}')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<FRAME>/{frame}/g" scan_job.sh > scan_job_mod.sh; rm scan_job.sh; mv scan_job_mod.sh scan_job.sh')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<FRAME>/{frame}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<HYDROGEN>/{H}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<OHATOM>/{O_OH1_atom}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<CXATOM>/{C_atom}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<HXATOM>/{H_atom}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan')}; sed "s/<FRAME>/{frame}/g" scan_plot.py > scan_plot_mod.py; rm scan_plot.py; mv scan_plot_mod.py scan_plot.py')
        os.mkdir(os.path.join(calculations, str(frame), '3_scan', 'structures'))
        execute(f'cp multiframe.py {os.path.join(calculations, str(frame), '3_scan', 'structures')}')
        execute(f'cp scan_view.py {os.path.join(calculations, str(frame), '3_scan', 'structures')}')
        execute(f'cd {os.path.join(calculations, str(frame), '3_scan', 'structures')}; sed "s:<PATH>:{os.path.abspath(os.path.join(calculations, str(frame), '3_scan', 'structures', 'multiframe_view.py'))}:g" scan_view.py > scan_view_mod.py; rm scan_view.py; mv scan_view_mod.py scan_view.py')

        # 4_TS_opt
        os.mkdir(os.path.join(calculations, str(frame), '4_TS_opt'))
        execute(f'cp {os.path.join(frame_source, str(frame), '*.prmtop')} {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cp {os.path.join(frame_source, str(frame), '*.inpcrd')} {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cp {os.path.join(frame_source, str(frame), '*.act_list')} {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cp {qm_list_name} {os.path.join(calculations, str(frame), '4_TS_opt', 'qm_list')}')
        execute(f'cp ts_info.txt {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cp ts_job.sh {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cp ts.chm {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cp opt_plot.py {os.path.join(calculations, str(frame), '4_TS_opt')}')
        execute(f'cd {os.path.join(calculations, str(frame), '4_TS_opt')}; mv opt_plot.py ts_plot.py')

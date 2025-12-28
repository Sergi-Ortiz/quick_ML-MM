#!/Users/sergiortizropero/miniconda3/envs/md_traj/bin/python
# Sergi Ortiz @ UAB
# 26 July 2025
# Set up QM/MM calculation via ChemShell @ B3LYP/6-31G(d)/AMBER(ff19SB)/TIP3P

# example usage:
# python setup_qmmm.py -H 10 -m 5-HETE



# 
#   IMPORT SECTION
#
import os
import shutil
import subprocess
import argparse


# macOS
bash = '/opt/homebrew/bin/bash'
# linux
#bash = '/bin/bash'

def execute(command_str):
    subprocess.run(command_str, shell=True, executable=bash)






# 
#   CLI ARGUMENTS
# 
parser = argparse.ArgumentParser(description='Set up QMMM calculations')
parser.add_argument(
    '-m', '--molecule', 
    type=str,
    required=True,
    help='Ligand used in the MD trajectory. Options: \'AA\', \'5-HETE\', \'5-HpETE\'.',
)
parser.add_argument(
    '-H', '--hydrogen', 
    type=int,
    required=True,
    help='Hydrogen selection for which to perform the QM/MM scan. Options: 10, 13.',
)
# 
#   QM/MM SETUP OPTIONS
#
args = parser.parse_args()
lig_ID = args.molecule
H = args.hydrogen

#
#   QM/MM INFORMATION
#
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
#   PATH DEFINITIONS
#
qm_list_name = f'qm_list_{lig_ID[2:]}'
FRAME_LIST = [x for x in os.listdir('../') if os.path.isdir(os.path.abspath(os.path.join('../', x))) and (x != 'qmmm_scripts' and x != 'results') ]

for frame in FRAME_LIST:
    BASE_DIR = os.path.join(f'../{frame}')
    SRC_DIR = os.path.join(BASE_DIR, 'src')
    QM_DIR = os.path.join(BASE_DIR, 'qm_region')
    OPT_DIR = os.path.join(BASE_DIR, 'opt')
    SCAN_DIR = os.path.join(BASE_DIR, 'scan')
    TS_DIR = os.path.join(BASE_DIR, 'ts')

    # QM_DIR generation
    if not os.path.exists(QM_DIR):
        os.mkdir(QM_DIR)
        execute(f'cp {qm_list_name} {QM_DIR}/qm_list')
        execute(f'cp {SRC_DIR}/*.pdb {QM_DIR}')
        execute(f'cp atoms_list_to_pdb.py {QM_DIR}; cd {QM_DIR}; python atoms_list_to_pdb.py \'qm_list\' {frame}_solv_cropped.pdb qm_list')
        
        # chimera script to quickly visualize the QM-region
        execute(f'cp qm_region_view.py {QM_DIR}; cd {QM_DIR}; ')
        execute(f'sed "s:<PATH>:{QM_DIR}/view.py:g" qm_region_view.py > qm_region_view_mod.py; rm qm_region_view.py; mv qm_region_view_mod.py qm_region_view.py')
    else:
        print(f'skipping qm_region generation, {QM_DIR} already exists.\n\n')


    # OPT_DIR generation
    if not os.path.exists(OPT_DIR):
        os.mkdir(OPT_DIR)

        # copy .prmtop .incprd .pdb and .act_list
        execute(f'cp {SRC_DIR}/* {OPT_DIR}')

        # copy qm_list
        execute(f'cp {QM_DIR}/qm_list {OPT_DIR}/qm_list')

        # copy opt scripts and modify them accordingly
        execute(f'cp opt_job.sh {OPT_DIR}; cp opt.chm {OPT_DIR}; cp plot.py {OPT_DIR}')
        execute(f'cd {OPT_DIR}; sed "s/<FRAME>/{frame}/g" opt_job.sh > opt_job_mod.sh; rm opt_job.sh; mv opt_job_mod.sh opt_job.sh')
        execute(f'cd {OPT_DIR}; sed "s/<FRAME>/{frame}/g" opt.chm > opt_mod.chm; rm opt.chm; mv opt_mod.chm opt.chm')
    else:
        print(f'skipping opt generation, {OPT_DIR} already exists.\n\n')


    # SCAN_DIR generation
    if not os.path.exists(SCAN_DIR):
        os.mkdir(SCAN_DIR)

        SCAN_DIR_STUCTURES = os.path.join(SCAN_DIR, 'structures')
        os.mkdir(SCAN_DIR_STUCTURES)

        # copy .prmtop .incprd and .act_list
        execute(f'cp {SRC_DIR}/* {SCAN_DIR}; cd {SCAN_DIR}; rm *.pdb')
        
        # copy qm_list
        execute(f'cp {QM_DIR}/qm_list {SCAN_DIR}/qm_list')

        # manual task
        # copy the resulting optimization opt.pdb to SCAN_DIR

        # copy scan scripts and modify them accordingly
        execute(f'cp scan_job.sh {SCAN_DIR}; cp scan_plot.py {SCAN_DIR}; cp scan.chm {SCAN_DIR}')
        execute(f'cd {SCAN_DIR}; sed "s/<FRAME>/{frame}/g" scan_job.sh > scan_job_mod.sh; rm scan_job.sh; mv scan_job_mod.sh scan_job.sh')
        execute(f'cd {SCAN_DIR}; sed "s/<FRAME>/{frame}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {SCAN_DIR}; sed "s/<HYDROGEN>/{H}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {SCAN_DIR}; sed "s/<OHATOM>/{O_OH1_atom}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {SCAN_DIR}; sed "s/<CXATOM>/{C_atom}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {SCAN_DIR}; sed "s/<HXATOM>/{H_atom}/g" scan.chm > scan_mod.chm; rm scan.chm; mv scan_mod.chm scan.chm')
        execute(f'cd {SCAN_DIR}; sed "s/<FRAME>/{frame}/g" scan_plot.py > scan_plot_mod.py; rm scan_plot.py; mv scan_plot_mod.py scan_plot.py')
        
        # add multiframe.py and scan_view.py to structures (where scan structures will be saved)
        execute(f'cp multiframe.py {SCAN_DIR_STUCTURES}')
        execute(f'cp scan_view.py {SCAN_DIR_STUCTURES}')
        execute(f'cd {SCAN_DIR_STUCTURES}; sed "s:<PATH>:{SCAN_DIR_STUCTURES}/multiframe_view.py:g" scan_view.py > scan_view_mod.py; rm scan_view.py; mv scan_view_mod.py scan_view.py')
    else:
        print(f'skipping scan generation, {SCAN_DIR} already exists.\n\n')


    # TS_SCAN generation
    if not os.path.exists(TS_DIR):
        os.mkdir(TS_DIR)

        # copy .inpcrd .prmtop and .act_list
        execute(f'cp {SRC_DIR}/* {TS_DIR}; cd {TS_DIR}; rm *.pdb')
        execute(f'cp {QM_DIR}/qm_list {TS_DIR}/qm_list')

        # copy TS scripts and modify them accordingly
        execute(f'cp ts_job.sh {TS_DIR}; cp ts.chm {TS_DIR}; cp plot.py {TS_DIR}; cp ts_info.txt {TS_DIR}')

        # manual task
        # read ts_info.txt, similarly to the optimization
    else:
        print(f'skipping ts generation, {TS_DIR} already exists.\n\n')
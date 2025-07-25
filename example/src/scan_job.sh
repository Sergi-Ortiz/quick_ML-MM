#!/bin/bash
#SBATCH -J <FRAME>_scan_qmmm
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.log
#SBATCH --mail-user=sergi.ortizr@autonoma.cat
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=30-00:00


# module load
ml Turbomole
ml Chemshell

# node computing options
export OMP_NUM_THREADS=16   
export PARNODES=16          
export PARA_ARCH=SMP
ulimit -s unlimited


# source dir    $SLURM_SUBMIT_DIR
# working dir   $TMP_DIR

# file management (input)
echo "<FRAME>    setting up input files"

cp $SLURM_SUBMIT_DIR/*<FRAME>* $TMP_DIR             # prmtop, inpcrd, c, pdb, act_list
cp $SLURM_SUBMIT_DIR/qm_list $TMP_DIR           # qm_list
cp $SLURM_SUBMIT_DIR/scan.chm $TMP_DIR           # ChemShell job
ls -l $TMP_DIR > $SLURM_SUBMIT_DIR/input_files.txt


# scan
echo "<FRAME>    starting QM/MM scan @ `date`"
cd $TMP_DIR/
srun chemsh scan.chm > $SLURM_SUBMIT_DIR/scan.out
echo "<FRAME>    finished QM/MM scan @ `date`"


# file management (output)
echo "<FRAME>    retrieving input files"
mkdir $SLURM_SUBMIT_DIR/output
cp $TMP_DIR/*.pdb $SLURM_SUBMIT_DIR/output
cp $TMP_DIR/*.c   $SLURM_SUBMIT_DIR/output
cp $TMP_DIR/*.plt $SLURM_SUBMIT_DIR

# file management (raw output)
echo "<FRAME>    retrieving raw files"
mkdir $SLURM_SUBMIT_DIR/output_raw
cp $TMP_DIR/* $SLURM_SUBMIT_DIR/output_raw

ls -l $TMP_DIR > $SLURM_SUBMIT_DIR/output/output_files.txt

echo "<FRAME>    finished job @ `date`"
echo ""
#!/bin/bash

#SBATCH -J 5-HETE_343
#SBATCH -e h15LOX_5-HETE_343_12_mod_norelaxed_4amber_solv.%j.err
#SBATCH -o h15LOX_5-HETE_343_12_mod_norelaxed_4amber_solv.%j.out
#SBATCH -p gorn1
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -t 30-00:00

# load Amber module
ml Amber

# copy all MD inputs to TEMP_DIR
cp -r * $TMP_DIR
cd $TMP_DIR


# set VARS
## general
prmtop=h15LOX_5-HETE_343_12_mod_norelaxed_4amber_solv.prmtop
inpcrd=h15LOX_5-HETE_343_12_mod_norelaxed_4amber_solv.inpcrd

# preproduction
declare -i preprod=1

# production
prefix='prod'
declare -i copy_preprod=1
declare -i prod=1
declare -i cnt=1
declare -i cntmax=25

if [ $preprod -eq 1 ];
then
    cd preprod

    # 1_min.in
    echo -e "step 1\t Starting first minimization step\t @ `date`"
    pmemd.cuda -O -i 1_min.in\
                -o out/1_min.out -p ../$prmtop -c ../$inpcrd -r out/1_min.rst\
                -inf out/1_min.info -ref ../$inpcrd -x out/1_min.nc
    echo -e "step 1\t Finished first minimization step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/
 
 
    # 2_heat.in
    echo -e "step 2\t Starting first heating step\t @ `date`"
    pmemd.cuda -O -i 2_heat.in\
                -o out/2_heat.out -p ../$prmtop -c out/1_min.rst -r out/2_heat.rst\
                -inf out/2_heat.info -ref out/1_min.rst -x out/2_heat.nc
    echo -e "step 2\t Finished first heating step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 3_npt.in
    echo -e "step 3\t Starting first NPT step\t @ `date`"
    pmemd.cuda -O -i 3_npt.in\
                -o out/3_npt.out -p ../$prmtop -c out/2_heat.rst -r out/3_npt.rst\
                -inf out/3_npt.info -ref out/2_heat.rst -x out/3_npt.nc
        echo -e "step 3\t Finished first NPT step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 4_npt.in
    echo -e "step 4\t Starting second NPT step\t @ `date`"
    pmemd.cuda -O -i 4_npt.in\
                -o out/4_npt.out -p ../$prmtop -c out/3_npt.rst -r out/4_npt.rst\
                -inf out/4_npt.info -ref out/3_npt.rst -x out/4_npt.nc
    echo -e "step 4\t Finished second NPT step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 5_min.in
    echo -e "step 5\t Starting second minimization step\t @ `date`"
    pmemd.cuda -O -i 5_min.in\
                -o out/5_min.out -p ../$prmtop -c out/4_npt.rst -r out/5_min.rst\
                -inf out/5_min.info -ref out/4_npt.rst -x out/5_min.nc
    echo -e "step 5\t Finished second minimization step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 6_npt.in
    echo -e "step 6\t Starting third NPT step\t @ `date`"
    pmemd.cuda -O -i 6_npt.in\
                -o out/6_npt.out -p ../$prmtop -c out/5_min.rst -r out/6_npt.rst\
                -inf out/6_npt.info -ref out/5_min.rst -x out/6_npt.nc
    echo -e "step 6\t Finished third NPT step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 7_npt.in
    echo -e "step 7\t Starting fourth NPT step\t @ `date`"
    pmemd.cuda -O -i 7_npt.in\
                -o out/7_npt.out -p ../$prmtop -c out/6_npt.rst -r out/7_npt.rst\
                -inf out/7_npt.info -ref out/6_npt.rst -x out/7_npt.nc
    echo -e "step 7\t Finished fourth NPT step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 8_npt.in
    echo -e "step 8\t Starting fifth NPT step\t @ `date`"
    pmemd.cuda -O -i 8_npt.in\
                -o out/8_npt.out -p ../$prmtop -c out/7_npt.rst -r out/8_npt.rst\
                -inf out/8_npt.info -ref out/7_npt.rst -x out/8_npt.nc
    echo -e "step 8\t Finished fifth NPT step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 9_npt.in
    echo -e "step 9\t Starting sixth NPT step\t @ `date`"
    pmemd.cuda -O -i 9_npt.in\
                -o out/9_npt.out -p ../$prmtop -c out/8_npt.rst -r out/9_npt.rst\
                -inf out/9_npt.info -ref out/8_npt.rst -x out/9_npt.nc
    echo -e "step 9\t Finished sixth NPT step\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/


    # 10_nvt.in
    echo -e "step 10\t Starting NVT equilibration\t @ `date`"
    pmemd.cuda -O -i 10_nvt.in\
                -o out/10_nvt.out -p ../$prmtop -c out/9_npt.rst -r out/10_nvt.rst\
                -inf out/10_nvt.info -ref out/9_npt.rst -x out/10_nvt.nc
    echo -e "step 10\t Finished NVT equilibration\t @ `date`"
    cp -r out $SLURM_SUBMIT_DIR/preprod/

    cd ..
fi

if [ $prod -eq 1 ];
then
    cd prod
    if [ $copy_preprod -eq 1 ]
    then
        cp ../preprod/out/10_nvt.rst .
    fi

    while [ ${cnt} -le ${cntmax} ]
    do
        pcnt=$((${cnt} - 1))
        istep=${prefix}_${cnt}
        pstep=${prefix}_${pcnt}

        if [ ${cnt} -eq 1 ]
        then
            # resume from preprod
            echo -e "Starting NVT production\t @ `date`"
            pmemd.cuda -O -ref 10_nvt.rst -p ../$prmtop -c 10_nvt.rst  -i ${istep}.in -o out/${istep}.mdout -inf ${istep}.mdinf -r out/${istep}.rst -x out/${istep}.nc
            echo -e "Finished first 10ns of NVT production\t @ `date`"
            cp -r out $SLURM_SUBMIT_DIR/prod/

        else
            # continue prod
            echo -e "Finished 10ns of NVT production\t @ `date`"
            pmemd.cuda -O -ref out/${pstep}.rst -p ../$prmtop -c out/${pstep}.rst  -i ${prefix}.in -o out/${istep}.mdout -inf ${istep}.mdinf -r out/${istep}.rst -x out/${istep}.nc

            echo -e "Finished 10ns of NVT production\t @ `date`"
	    cp -r out $SLURM_SUBMIT_DIR/prod/
        fi

        cnt=$(($cnt + 1))
        echo "Step ${cnt} of ${cntmax} of production done"
    done
fi


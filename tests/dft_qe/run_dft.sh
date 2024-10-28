##!/bin/bash -l
##SBATCH -J "X_BSE_s"
##SBATCH -N 1
##SBATCH --ntasks-per-node=28
##SBATCH --cpus-per-task=1
##SBATCH --time=48:00:00
##SBATCH --hint=multithread
##SBATCH --qos=normal


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK


#module load toolchain/intel/2020b

#export PATH="$PATH:/home/users/mnalabothula/iris/q-e/bin"


MPI_CMD=mpirun
export SLURM_NTASKS=4

#MPI_CMD=srun
PYTHON=python3
PW=pw.x
PH=ph.x
P2Y=p2y
YAMBO=yambo
LELPHC=lelphc


NCPUS=$SLURM_NTASKS
NKpool=1

test_systems=(\
    Si_bulk_spin \
    Si_bulk \
    hBN_bulk \
    NiO_SOC_mag \
    MoS2_SOC_2D \
    MoS2_SOC_2D_no_soc1 \
    MoS2_SOC_2D_no_soc2 \
    MoS2_SOC_2D_no_soc3 \
    )


for test_str in ${test_systems[@]}
do
    cd "$test_str"
    rm -rf scf phonons nscf 
    mkdir scf phonons nscf 
    cp scf.in scf
    cd scf

    # scf
    $MPI_CMD -n $NCPUS $PW < scf.in | tee scf.out
    cd ../phonons
    cp ../ph.in .
    cp -r ../scf/*.save .
    cp -r ../scf/*.xml .

    ## phonons
    $MPI_CMD -n $NCPUS $PH  <ph.in | tee ph.out
    $MPI_CMD -n $NCPUS $LELPHC -pp --code=qe -F ph.in 
    rm *wfc*
    cd ../nscf
    cp ../nscf.in .
    cp -r ../scf/*.save .
    cp -r ../scf/*.xml .

    #symmetric
    $MPI_CMD -n $NCPUS $PW  < nscf.in | tee nscf.out
    rm *wfc*
    cd *.save
    $MPI_CMD -n 1 $P2Y
    $MPI_CMD -n 1 $YAMBO
    cd ../..
    cp -r nscf/*.save/SAVE .
    cp -r phonons/ph_save .
    rm -rf scf phonons nscf

    ## Convert to portable databases
    # 1) COnvert ph_save 
    $PYTHON ../../convert_data.py --to_npy=ph_save 
    # 2) convert SAVE to float
    $PYTHON ../../convert_data.py --to_float=SAVE 
    cd ..
done


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
PW=pw.x
PH=ph.x
P2Y=p2y
YAMBO=yambo
LELPHC=lelphc

NCPUS=$SLURM_NTASKS
NKpool=1

test_systems=(\
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
    cd phonons
    $MPI_CMD -n $NCPUS $LELPHC -pp --code=qe -F ph.in
    cd ..
    cd nscf
    cp ../../input.in .
    ln -sf *.save/SAVE .
    $MPI_CMD -n $NCPUS $LELPHC --code=qe -F input.in
    cd ..
    cd nscf_wo
    cp ../../input.in .
    ln -sf *.save/SAVE .
    $MPI_CMD -n $NCPUS $LELPHC --code=qe -F input.in
    cd ../..
done


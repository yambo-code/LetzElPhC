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
    rm -rf scf phonons nscf nscf_wo
    mkdir scf phonons nscf nscf_wo
    cp scf.in scf
    cd scf

    # scf
    $MPI_CMD -n $NCPUS $PW < scf.in | tee scf.out
    cd ../phonons
    cp ../ph.in .
    cp -r ../scf/*.save .
    cp -r ../scf/*.xml .

    ## phonons
    $MPI_CMD -n $NCPUS $PH  < ph.in | tee ph.out
    rm *wfc*
    cd ../nscf
    cp ../nscf.in .
    cp -r ../scf/*.save .
    cp -r ../scf/*.xml .

    #symmetric
    $MPI_CMD -n $NCPUS $PW  < nscf.in | tee nscf.out
    rm *wfc*
    cd *.save
    cp ../../dipoles.in .
    $MPI_CMD -n 1 $P2Y
    $MPI_CMD -n 1 $YAMBO
    $MPI_CMD -n $NCPUS $YAMBO -F dipoles.in # dipoles
    ### elph for symmetric 
    cd ..
    cp ../elph.in .
    cp -r ../phonons/_ph* ../phonons/*.dyn* .
    #$MPI_CMD -n $NCPUS $PH  < elph.in | tee elph.out
    #cd elph_dir
    #cp ../../../convert_scripts/qe2yambo .
    #echo "s.dbph_000001" | ./qe2yambo
    #cd ..

    cd ../nscf_wo
    cp ../nscf_wo.in .
    cp -r ../scf/*.save .
    cp -r ../scf/*.xml .

    ##no_sym
    $MPI_CMD -n $NCPUS $PW < nscf_wo.in | tee nscf.out
    rm *wfc*
    cd *.save
    cp ../../dipoles.in .
    $MPI_CMD -n 1 $P2Y
    $MPI_CMD -n 1 $YAMBO
    $MPI_CMD -n $NCPUS $YAMBO -F dipoles.in
    cd ..
    cp ../elph.in .
    cp -r ../phonons/_ph* ../phonons/*.dyn* .
    #$MPI_CMD -n $NCPUS $PH  < elph.in | tee elph.out
    #cd elph_dir
    #cp ../../../convert_scripts/qe2yambo .
    #echo "s.dbph_000001" | ./qe2yambo
    #cd ..
    cd ../..
done


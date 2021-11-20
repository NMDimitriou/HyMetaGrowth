#!/bin/bash
#SBATCH --account=
#SBATCH --job-name=CAL_KSC_FW_SET_2_RUN_1_IC_EXP
#SBATCH --output=%x-%j.out
#SBATCH --time=0-24:00
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --nodes=16
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --mem=10G

module load nixpkgs/16.09 gcc/7.3.0  gsl/2.5 cuda/10.1 openmpi/3.1.2
scontrol show hostnames > nodelist-$SLURM_JOB_ID
export PATH=$HOME/usr/torc/bin:$PATH

mpirun --map-by ppr:1:node ./main_parent start ../../IC/series2/Control_s2_FW Control_s2_FW
cd ../../

sleep 5

mpirun -hostfile nodelist-$SLURM_JOB_ID ./sample
#mpirun -np 64 ./sample

sleep 5

cd model/parent/
mpirun --map-by ppr:1:node ./main_parent stop ../../IC/series2/Control_s2_FW Control_s2_FW
cd ../../

mkdir $SLURM_JOB_NAME
mv *.txt *.out $SLURM_JOB_NAME

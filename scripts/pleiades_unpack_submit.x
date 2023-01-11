#!/bin/bash
#PBS -W group_list=s2317
#PBS -l select=7:ncpus=18:model=bro
###PBS -q 
#PBS -q devel
#PBS -l walltime=00:20:00
#PBS -j oe
#PBS -W umask=33
#PBS -m bea
#PBS -r n
####
#PBS -N CCSHIGHRES_unpack 
#PBS -o CCSHIGHRES_unpack.out 
#PBS -e CCSHIGHRES_unpack.err 
#PBS -V
####

module purge
module load comp-intel/2018.3.222 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
####module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
####module load comp-intel/2018.3.222  mpi-hpe/mpt.2.21 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
module list


ASSIM_RUNDIR=/nobackup/averdy/CCS/run_ad/
cd $ASSIM_RUNDIR
echo $PWD

cp mitgcmuv_unpack mitgcmuv
cp data.ctrl_pack data.ctrl

#set up optimization iteration number

mpiexec -np 126 /u/scicon/tools/bin/mbind.x ./mitgcmuv > run_unpack.out 


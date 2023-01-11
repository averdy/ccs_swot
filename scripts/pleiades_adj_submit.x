#!/bin/bash
#PBS -W group_list=s2317
#PBS -l select=7:ncpus=18:model=bro
###PBS -q devel
#PBS -q long
#PBS -l walltime=30:00:00
#PBS -j oe
#PBS -W umask=33
#PBS -m bea
#PBS -r n
####
#PBS -N CCSHIGHRES_adj 
#PBS -o CCSHIGHRES_adj.out 
#PBS -e CCSHIGHRES_adj.err 
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

cp data_ad data
sed -i "s|TEMPDIR|$TMPDIR|g" data
cp data.pkg_ad data.pkg
cp data.ctrl_ad data.ctrl
cp data.ecco_ad data.ecco

#set up optimization iteration number

mpiexec -np 126 /u/scicon/tools/bin/mbind.x ./mitgcmuv_ad > run_for_adj.out 


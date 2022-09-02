CCS state estimation on Pleiades
--------------------------------------
Ariane V. (09/01/22)
Based on notes from Ganesh G.
and code from Alex A.


Summary:
The high-resolution CCS regional state estimate runs on the NASA Pleiades computer. The run is for October 1-22, 2019 (22 days) and constraints are altimetry, satellite SST, and Argo profiles. The model is forced with ERA5, with initial and open boundary conditions from HYCOM reanalysis. Two iterations are run (forward+adjoint, packing, optim, unpacking) and the cost descends by ***%. To reproduce this, follow the steps below.


-----------------
# Model code

1) obtain MITgcm checkpoint 68i <br />
% git clone https://github.com/MITgcm/MITgcm.git <br />
% cd MITgcm <br />
% git checkout checkpoint68i


2) Obtain CCS-specific code and inputs <br />
% git clone https://github.com/averdy/ccs_swot <br />


Notes on the code:
- The size of the model grid is 774x966, with 126 tiles of size 86x69. Thus it requires 126 processes to run. With Alex's original 36-tile configuration, the memory use is too high (nodes crash). <br />
- The time step is 300 seconds (5 minutes). It takes *** minutes to run 22 days fwd+adjoint (45 minutes fwd). <br />
- It is set up for a 22-day run. For longer runs, need to adjust tamc.h  <br />
- Still tuning CPP_OPTIONS, etc <br />


-----------------
# Compiling 

1) Loading modules <br />
Modify ~/.bashrc to look like this: <br />

------------------  <br />
.bashrc <br />

if [ -f /etc/bashrc ]; then <br />
        . /etc/bashrc <br />
fi <br />

export PATH=/home4/averdy/STAF:${PATH} <br />

Load modules module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt  <br />netcdf/4.4.1.1_mpt <br />
------------------  <br />


2) Set up TAF

Copy  taf.pub  keys (from .ssh/ on mist.ucsd.edu) to /home4/averdy/.ssh/ 

testing TAF: <br />
% staf -test <br />


3) Create blas libraries

Compile ~/MITgcm/lsopt

Modify Makescript, replacing  <br />
FC              = f77 <br />
FFLAGS          = -fconvert=big-endian -fimplicit-none <br />
with <br />
FC              = ifort <br />
FFLAGS          = -mcmodel=large -shared-intel -fp-model precise -132 -r8 -i4 -W0 -WB -CB -fpe0 -traceback -convert big_endian -assume byterecl <br />

and then use   <br />
% make all <br />


4) Compile the code 

Generate executable mitgcmuv_ad: <br />
cd /home4/averdy/MITgcm/assim/CCS/build_ad/ <br />
./makescript_adj.sio.pleiades  <br />

The makescript calls options that are in <br />
/home4/averdy/MITgcm/assim/pleiades_build_options/linux_amd64_ifort+mpi_ice_nas  <br />
and includes the necessary MPI headers <br />

Generate pack/unpack executables: <br />
% cd /home4/averdy/MITgcm/assim/CCS/build_pack/ <br />
% ./makescript_fwd.sio.pleiades  <br />
% cp mitgcmuv mitgcmuv_pack <br />
% cd /home4/averdy/MITgcm/assim/CCS/build_unpack/ <br />
% ./makescript_fwd.sio.pleiades  <br />
% cp mitgcmuv mitgcmuv_unpack <br />


5) Compile line-search algorithm <br />

Generate optim.x in /home4/averdy/MITgcm/optim/ <br />

Modify Makescript, replacing  <br />
INCLUDEDIRS     = -I.                           \ <br />
                  -I../verification/tutorial_global_oce_optim/build/ <br />
with  <br />
INCLUDEDIRS     = -I.                           \ <br />
                  -I../assim/CCS/build_pack/ <br />
and <br />
       -DMAX_INDEPEND=1000000          \ <br />
with  <br />
       -DMAX_INDEPEND=305584550        \ <br />
and <br />
LIBS            = -llsopt_ecco                 \ <br />
                  -lblas1 <br />
with <br />
LIBS            = -llsopt_ecco                  \ <br />
                  -mkl                          \ <br />
                  -lpthread                     \ <br />
                  -lmpi <br />
and  <br />
FC              = f77 <br />
FFLAGS          = -fconvert=big-endian -fimplicit-none <br />
with  <br />
FC              = ifort <br />
FFLAGS          =  -mcmodel=large -shared-intel -fp-model precise -132 -r8 -i4 -W0 -WB -CB -fpe0 -traceback -convert big_endian -assume byterecl <br />

and then use <br />
% make depend <br />
% make <br />


-----------------
# Model inputs

1) grid: <br />
- Bathymetry is from Alex: TFO_2km_bathy.bin <br />
- Vertical resolution is Ariane's 100 levels: delRFile_100_5100m.bin <br />
- Run for a few time steps to generate grid files (XC.data, YC.data, etc) <br />

2) constraints: <br />
- Obtain Argo profiles for 2019 (processed by Sharon E.) <br />
- Make SST constraint for 2019 (make_OISST_ccs.m) <br />
- Make SSH constraint for 2019 (Krads2grd_CCS.E.PER_YEAR / rads_QC_peryear.m) <br />
- Make geoid constraint (mdt_products_regrid.m) <br />

3) weights: <br />
- Make uncertainties for ICS and SST (make_mapped_weights_from_Argo_error.m), <br />
atm state (ERA_make_weights.m), <br />
SSH (make_error_ssh.m; 3 or 6 cm uniform), <br />
geoid (10 cm uniform) <br />


4) forcing: <br />
- Obtain ERA5 for 2019 <br />
- Make ICs from HYCOM - nov 1, 2019 (make_ics_hycom.m) <br />
- Make OBCs from HYCOM - jan-dec 2019 (make_obcs_hycom.m) <br />
- Make runoff from CORE climatology (make_runoff_core.m) <br />
- *** tides <br />


5) MITgcm data files: <br />
- Take Alex's files, from LLC4320, and change: <br />
- ...


-----------------
# Running

Create run directory <br />
% mkdir /nobackup/averdy/CCS/run_ad <br />



-----------------


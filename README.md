CCS state estimation on Pleiades
--------------------------------------
Ariane V. (09/01/22)
Based on notes from Ganesh G.
and code from Alex A.


Summary:
The high-resolution CCS regional state estimate runs on the NASA Pleiades computer. The run is for October 1-22, 2019 (22 days) and constraints are satellite SST. Currently, profiles, geoid and altimetry constraints are not working. The model is forced with ERA5, with initial and open boundary conditions from HYCOM reanalysis. Two iterations are run (adjoint, packing, optim, unpacking) and the cost descends xxx%. To reproduce this, follow the steps below.


*** Model code

1) obtain MITgcm checkpoint 68i <br />
% git clone https://github.com/MITgcm/MITgcm.git <br />
% cd MITgcm <br />
% git checkout checkpoint68i


2) Obtain CCS-specific code and inputs <br />
% git clone https://github.com/... 


Notes on the code:
- The size of the model grid is 774x966, with 126 tiles of size 86x69. Thus is requires 126 processes to run. I also tried a 36-tile configuration, but the memory use was too high (nodes crash). <br />
- The time step is 300 seconds (5 minutes). It takes * minutes to run 22 days fwd+adjoint (45 minutes fwd). <br />
- It is set up for a 22-day run. For longer runs, need to adjust tamc.h  <br />
- Still tuning CPP_OPTIONS, etc <br />


*** Compiling 

CCS state estimation on Pleiades
--------------------------------------
Ariane (08/25/22)
Based on notes from Ganesh
and code from Alex


Summary:
The high-resolution CCS regional state estimate runs on the NASA Pleiades computer. The run is for October 1-22, 2019 (22 days) and constraints are satellite SST. Currently, profiles, geoid and altimetry constraints are not working. The model is forced with ERA5, with initial and open boundary conditions from HYCOM reanalysis. Two iterations are run (adjoint, packing, optim, unpacking) and the cost descends xxx%. To reproduce this, follow the steps below.


*** MITgcm code

1) obtain MITgcm checkpoint 68i <br />
git clone https://github.com/MITgcm/MITgcm.git <br />
cd MITgcm <br />
git checkout checkpoint68i


2) Obtain CCS-specific code and inputs <br />

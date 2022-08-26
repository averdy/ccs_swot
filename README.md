CCS STATE ESTIMATION ON PLEIADES
--------------------------------------
Ariane (08/25/22)
Based on notes from Ganesh
and code from Alex


Summary:
The high-resolution CCS regional state estimate runs on the NASA Pleiades computer. The run is for October 1-22, 2019 (22 days) and constraints is satellite SST. Currently, profiles, geoid and altimetry constraints are not working. The model is initialized with HYCOM state and forced with ERA5, with HYCOM open boundary conditions. Two iterations are run (adjoint, packing, optim, unpacking) and the cost descends xxx%. Details on how to produce the state estimate are documented below.


*** MITgcm code

1) obtain MITgcm checkpoint 68i <br />
git clone https://github.com/MITgcm/MITgcm.git <br />
cd MITgcm <br />
git checkout checkpoint68i


2) Obtain CCS-specific code and inputs <br />

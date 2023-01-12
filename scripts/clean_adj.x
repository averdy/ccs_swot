
rm /nobackup/averdy/CCS/run_ad/mitgcm*
rm /nobackup/averdy/CCS/run_ad/cost*
rm /nobackup/averdy/CCS/run_ad/diag*
rm /nobackup/averdy/CCS/run_ad/tape*
rm /nobackup/averdy/CCS/run_ad/PROF/*
rm /nobackup/averdy/CCS/run_ad/STD*
rm /nobackup/averdy/CCS/run_ad/*meta
rm /nobackup/averdy/CCS/run_ad/*data
rm /nobackup/averdy/CCS/run_ad/scratch*
rm /nobackup/averdy/CCS/run_ad/*nc
rm /nobackup/averdy/CCS/run_ad/ecco* 
rm /nobackup/averdy/CCS/run_ad/*out
rm /nobackup/averdy/CCS/run_ad/OPW*

cp build_ad/mitgcmuv_ad /nobackup/averdy/CCS/run_ad
cp build_pack/mitgcmuv /nobackup/averdy/CCS/run_ad/mitgcmuv_pack
cp build_unpack/mitgcmuv /nobackup/averdy/CCS/run_ad/mitgcmuv_unpack
cp ../../optim/optim.x /nobackup/averdy/CCS/run_ad/

cp input_ad/* /nobackup/averdy/CCS/run_ad

ln -sf /nobackup/hzhang1/forcing/era5/ERA5* /nobackup/averdy/CCS/run_ad/
ln -sf /nobackup/averdy/CCS/SWOT_input/* /nobackup/averdy/CCS/run_ad/

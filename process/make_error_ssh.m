% constant error, 0.02 m

cd /data/SO6/CCS/
load /data/SO6/CCS/grid/grid XC YC Depth
err=0*XC+0.02;
err(Depth==0)=0;
fid = fopen('ssh_error_ccs_2cm.bin','w','b');
fwrite(fid,err,'single');
fclose(fid);



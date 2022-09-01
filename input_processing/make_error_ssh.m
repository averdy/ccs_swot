% constant error, 0.03 m or 0.06 m

cd /data/SO6/CCS/
load /data/SO6/CCS/grid/grid XC YC Depth
err=0*XC+0.03;
err(Depth==0)=0;
fid = fopen('ssh_error_ccs_3cm.bin','w','b');
fwrite(fid,err,'single');
fclose(fid);

err=err*2;
fid = fopen('ssh_error_ccs_6cm.bin','w','b');                                                                                                                      
fwrite(fid,err,'single');                                                                                                                                          
fclose(fid);                                                                                                                                                       


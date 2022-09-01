
outpath = '/data/SO6/CCS/';


% Rio 13 

cd /project_shared/cnes_cls2013

lat = ncread('mdt_cnes_cls2013_global.nc','lat');
lon = ncread('mdt_cnes_cls2013_global.nc','lon');
mdt = ncread('mdt_cnes_cls2013_global.nc','mdt');

load /data/SO6/CCS/grid/grid XC YC Depth

mdt_interp = interp2(lat,lon,mdt,YC,XC);
mdt_interp(Depth==0)=-9999;
mdt_interp(isnan(mdt_interp))=-9999;

% cm
mdt_interp = mdt_interp*100;

fid=fopen([outpath '/mdt_cnes_cls2013_cm_ccs.bin'],'w','b');
fwrite(fid,mdt_interp,'single');
fclose(fid);



% DTU 17
% DTU 19

cd /data/averdy/datasets/geoid

% read ascii file
%A = dlmread('dtu17mdt2.grd');
A = dlmread('dtu19mdt.grd');

% first line is grid information
lat0=A(1,1); lat1=A(1,2); dlat=A(1,5);
lon0=A(1,3); lon1=A(1,4); dlon=A(1,6);
lat = lat0:dlat:lat1;
lon = lon0:dlon:lon1;

mdt = A(2:end,:)';

% not sure why, but mdt read goes from 90N to 80S
lat = flip(lat);
% figure;pcolor(lon,lat,mdt');shading flat;colorbar

mdt_interp = interp2(lat,lon,mdt,YC,XC);
mdt_interp(Depth==0)=-9999;
mdt_interp(isnan(mdt_interp))=-9999;

% cm
mdt_interp = mdt_interp*100;

fid=fopen([outpath '/dtu19mdt2_cm_ccs.bin'],'w','b');
fwrite(fid,mdt_interp,'single');
fclose(fid);



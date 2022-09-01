% read HYCOM expt 93.0, 2018-2020
% to make obcs for CCS domain

clear
addpath ~/scripts_m

datestr = '2019';

% open data file
% Jan-1-2018 to Feb-19-2020
OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0';
ncid = netcdf.open(OpenDAP_URL,'NOWRITE');
prec = 'double';

% Variable id in data set
nc_depth = 4;
nc_date = 2;
nc_lat = 0;
nc_lon = 1;
nc_ssh = 5;
nc_salinity = 8; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,8)
nc_temperature = 6; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,6)
nc_u = 10; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,10)
nc_v = 12; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,12)

% time
% 3-hourly output, 5-day sampling
dt = 5; dthycom=dt*(24/3);
hycom_date = netcdf.getVar(ncid,nc_date);
% units are hours since 2000-1-1
hycom_date = datenum(2000,1,1)+hycom_date/24;
% my OBCS starting in 2019
obcs_date = datenum(2019,1,1):dt:(datenum(2019,12,31)+dt);
% first record to download 
rec1 = hycom_date(find(hycom_date>=obcs_date(1),1));
tmin = find(hycom_date==rec1);
rec2 = hycom_date(find(hycom_date>=obcs_date(end),1));
tmax = find(hycom_date==rec2);
date = hycom_date(tmin):dt:(hycom_date(tmax));
nt = length(date);
tmin = tmin-1; % for data download
tmax = tmax-1; 
%this doesn't work because missing days
%hycom_date = netcdf.getVar(ncid,nc_date,tmin,nt,dthycom,prec);
% !!!!!! temporary
% because loading nt records is not working
%nt = nt-1;
%hycom_date = hycom_date(1:end-1);

nthycom=length(hycom_date(tmin:tmax));
% hycom_date = netcdf.getVar(ncid,nc_date,tmin,nthycom,prec);


% There are missing days (!)
% read 1 record per day
ind = [];
days = (date(1)-2):date(end)+2;
for t=1:length(days)
% skip if doens't exist
 tmpind = find(hycom_date==days(t));
 if tmpind>0
  ind(t) = tmpind;
 else
  ind(t) = NaN;
  display(['date missing, t=' num2str(t)])
 end
end



 % 230E to 246.2E;
 xmin = 2875; xmax = 3079; nx = xmax-xmin+1;
 % 27N to 43.1N
 ymin = 1838; ymax = 2079; ny = ymax-ymin+1;
 zmin = 0; zmax = 39; nz = zmax-zmin+1;


hycom_depth = netcdf.getVar(ncid,nc_depth,prec);
hycom_depth = -1*hycom_depth;
% HYCOM GOES ONLY UPTO 5000 m, so adding one more layer
hycom_depth(end+1) = -6500;

hycom_lat = netcdf.getVar(ncid,nc_lat,ymin,ny,prec);


hycom_lon = netcdf.getVar(ncid,nc_lon,xmin,nx,prec);
%if hycom_lon < 0; hycom_lon = hycom_lon + 360; end



% read dynamical fields
params = {'V', 'U', 'T', 'S'};
ctl = {nc_v,nc_u,nc_temperature,nc_salinity};

for np = 1:length(ctl)

clear obc* 

nc_var = ctl{np};
param = params{np}


cnt=0;
for t=1:length(ind)
cnt=cnt+1
if ind(t)>0
if np==1
 obcw(:,:,cnt) = zeros(ny,nz);
else
 tmp = (netcdf.getVar(ncid,nc_var,[xmin,ymin,zmin,ind(t)],[1,ny,nz,1],prec)); 
 obcw(:,:,cnt) = squeeze(tmp(1,:,:));
end
if np==2
 obcs(:,:,cnt) = zeros(nx,nz);
 obcn(:,:,cnt) = zeros(nx,nz);
else
 tmp = (netcdf.getVar(ncid,nc_var,[xmin,ymin,zmin,ind(t)],[nx,1,nz,1],prec));
 obcs(:,:,cnt) = squeeze(tmp(:,1,:));
 tmp = (netcdf.getVar(ncid,nc_var,[xmin,ymax,zmin,ind(t)],[nx,1,nz,1],prec));
 obcn(:,:,cnt) = squeeze(tmp(:,1,:));
end
else
obcw(:,:,cnt) = NaN(ny,nz);
obcn(:,:,cnt) = NaN(nx,nz);
obcs(:,:,cnt) = NaN(nx,nz);
end
end

%obcw(:,end+1,:) = obcw(:,end,:);
%obcn(:,end+1,:) = obcn(:,end,:);
%obcs(:,end+1,:) = obcs(:,end,:);


% scaling
% from Ganesh
scale = [0.001, 0.001, 0.001, 0.001];
offset = [0, 0, 20, 20];


miss_value = -3e4;
obcs(obcs <= miss_value) = nan;
obcn(obcn <= miss_value) = nan;
obcw(obcw <= miss_value) = nan;

obcs=obcs*scale(np)+offset(np);
obcn=obcn*scale(np)+offset(np);
obcw=obcw*scale(np)+offset(np);



save(['/data/SO6/CCS/obcs/HYCOM_raw_daily/hycom_' param '_' datestr],'obcs','obcn','obcw');

end


% make OBCS by interpolating hycom data onto model grid
% (download hycom data first using get_hycom_ccs_obcs.m)
% - adjust transport so that it matches hycom transport after interpolation
% - balance transport by changing obcw so that net transport is 0
% - balance transport by changing obcw so that mean SSH varies like AVISO


datestr = '2019';

addpath ~/scripts_m

% read model grid

load /data/SO6/CCS/grid/grid XC YC RC DXG DYG DRF hFacS hFacW hFacC RAC
[nx_obc ny_obc nz_obc]=size(hFacC);


cd('/data/SO6/CCS/obcs/');


% for calculating hycom volume transport

% open data file
% Jan-1-2018 to Feb-2020
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

 % 230E to 246.2E;
 xmin = 2875-1; xmax = 3079-1; nx = xmax-xmin+1;
 % 27N to 43.1N
 ymin = 1838-1; ymax = 2079-1; ny = ymax-ymin+1;
 zmin = 0; zmax = 39; nz = zmax-zmin+1;

hycom_depth = netcdf.getVar(ncid,nc_depth,prec);
hycom_depth = -1*hycom_depth;
% HYCOM GOES ONLY UPTO 5000 m, so adding one more layer
hycom_depth(end+1) = -6500;

hycom_lat = netcdf.getVar(ncid,nc_lat,ymin,ny,prec);


hycom_lon = netcdf.getVar(ncid,nc_lon,xmin,nx,prec);
%if hycom_lon < 0; hycom_lon = hycom_lon + 360; end

% to calculate transport in hycom
dz = diff(hycom_depth(1:end-1)); dz = -[dz; dz(end)];
dx = diff(hycom_lon); dx = [dx; dx(end)]; dx = dx*111111; % convert to meters
dy = diff(hycom_lat); dy = [dy; dy(end)]; dy = dy*111111; % convert to meters
% area
for j=1:ny
for k=1:nz
dA_west(j,k) = dy(j)*dz(k);
end
end
for i=1:nx
for k=1:nz
dA_south(i,k) = dx(i)*cos(hycom_lat(1)*pi/180)*dz(k);
dA_north(i,k) = dx(i)*cos(hycom_lat(end)*pi/180)*dz(k);
end
end



% make obc every 5 days
dt_obc = 5;
my_obc_date = datenum(2019,1,1):dt_obc:datenum(2020,1,4);
nt_obc = length(my_obc_date);

% hycom dates (saved 1 record per day)
hycom_dates = datenum(2019,1,1):datenum(2020,1,4);


% loop parameters

params = {'V', 'U', 'T', 'S'};


for np = 1:length(params)

param = params{np}

tmpstr={'2019'};
obcn_all=[];
obcs_all=[];
obcw_all=[];
nt1=1;

% read all data for that parameter
for tt=1
load(['/data/SO6/CCS/obcs/HYCOM_raw_daily/hycom_' param '_' char(tmpstr{tt})]);
nt2=size(obcn,3);
% prep for regridding
obcn_all(:,:,nt1:nt1+nt2-1) = obcn; 
obcs_all(:,:,nt1:nt1+nt2-1) = obcs; 
obcw_all(:,:,nt1:nt1+nt2-1) = obcw; 
nt1=nt1+nt2;
end
nt=nt1-1; clear nt1 nt2
clear obcn obcs obcw


% repeat bottom layer
obcn_all(:,nz+1,:)=obcn_all(:,nz,:);
obcs_all(:,nz+1,:)=obcs_all(:,nz,:);
obcw_all(:,nz+1,:)=obcw_all(:,nz,:);


for t=1:nt_obc
 % read 2 days before and after to make 5-day average
 % except first time
 if t==1
  dates_to_read = my_obc_date(t):my_obc_date(t)+2;
 else
  dates_to_read = my_obc_date(t)-2:my_obc_date(t)+2;
 end

clear ind 
for d=1:length(dates_to_read)
 ind(d) = find(dates_to_read(d)==hycom_dates);
end

obcn(:,:,t)=nanmean(obcn_all(:,:,ind),3);
obcs(:,:,t)=nanmean(obcs_all(:,:,ind),3);
obcw(:,:,t)=nanmean(obcw_all(:,:,ind),3);



% calculate hycom volume transport

if strcmp(param,'V')

Tr_n_hycom(t) = nansum(nansum(obcn(:,1:nz,t).*dA_north));
Tr_s_hycom(t) = nansum(nansum(obcs(:,1:nz,t).*dA_south));

elseif strcmp(param,'U')

Tr_w_hycom(t) = nansum(nansum(obcw(:,1:nz,t).*dA_west));

end

end

clear obc*all




% prep for regridding

obcn_tmp=obcn;
obcs_tmp=obcs;
obcw_tmp=obcw;

for repeat=1:10
for k=1:nz
for i=2:nx-1
if isnan(obcs_tmp(i,k,1))
obcs_tmp(i,k,:) = nanmean(obcs_tmp(i-1:i+1,k,:));
end
if isnan(obcn_tmp(i,k,1))
obcn_tmp(i,k,:) = nanmean(obcn_tmp(i-1:i+1,k,:));
end
end
for i=2:ny-1
if isnan(obcw_tmp(i,k,1))
obcw_tmp(i,k,:) = nanmean(obcw_tmp(i-1:i+1,k,:));
end
end
end
end

for i=1:nx
for k=2:nz+1
if isnan(obcs_tmp(i,k,1))
obcs_tmp(i,k,:) = obcs_tmp(i,k-1,:);
end
if isnan(obcn_tmp(i,k,1))
obcn_tmp(i,k,:) = obcn_tmp(i,k-1,:);
end
end
end
for i=1:ny
for k=2:nz+1
if isnan(obcw_tmp(i,k,1))
obcw_tmp(i,k,:) = obcw_tmp(i,k-1,:);
end
end
end


if np==1 % V
maskS = squeeze(hFacS(:,2,:)); maskS(maskS>0)=1;
maskN = squeeze(hFacS(:,end,:)); maskN(maskN>0)=1;
maskW = squeeze(hFacS(1,:,:)); maskW(maskW>0)=1;
elseif np==2 % U
maskS = squeeze(hFacW(:,1,:)); maskS(maskS>0)=1;
maskN = squeeze(hFacW(:,end,:)); maskN(maskN>0)=1;
maskW = squeeze(hFacW(2,:,:)); maskW(maskW>0)=1;
else
maskS = squeeze(hFacC(:,1,:)); maskS(maskS>0)=1;
maskN = squeeze(hFacC(:,end,:)); maskN(maskN>0)=1;
maskW = squeeze(hFacC(1,:,:)); maskW(maskW>0)=1;
end



% regrid
tmp_s=0; tmp_n=0;tmp_w=0;
for t=1:nt_obc
tmp = interp2(hycom_lon,hycom_depth,obcs_tmp(:,:,t)',XC(:,1)',RC)';
tmp(isnan(tmp))=0;			
obcs_reg(:,:,t) = tmp.*maskS;
tmp_s=tmp_s+sum(sum(obcs_reg(:,:,t)==0&maskS>0));			
tmp = interp2(hycom_lon,hycom_depth,obcn_tmp(:,:,t)',XC(:,1)',RC)';
tmp(isnan(tmp))=0;			
obcn_reg(:,:,t) = tmp.*maskN;
tmp_n=tmp_n+sum(sum(obcn_reg(:,:,t)==0&maskN>0));
tmp = interp2(hycom_lat,hycom_depth,obcw_tmp(:,:,t)',YC(1,:),RC)';
tmp(isnan(tmp))=0;			
obcw_reg(:,:,t) = tmp.*maskW;
tmp_w=tmp_w+sum(sum(obcw_reg(:,:,t)==0&maskW>0));
end
tmp_s,tmp_n,tmp_w


if strcmp(param,'V')
for t=1:nt_obc

maskN = squeeze(hFacS(:,end,:));
maskS = squeeze(hFacS(:,2,:));
areaN = ((DXG(:,end)*DRF').*maskN);
areaS = ((DXG(:,2)*DRF').*maskS);
Tr_n_tmp = nansum(nansum(obcn_reg(:,:,t).*areaN));
Tr_s_tmp = nansum(nansum(obcs_reg(:,:,t).*areaS));
tmp = obcn_reg(:,:,t);
tmp(maskN>0) = tmp(maskN>0) - (Tr_n_tmp - Tr_n_hycom(t))/sum(sum(areaN));
obcn_corr(:,:,t) = tmp;
tmp = obcs_reg(:,:,t);
tmp(maskS>0) = tmp(maskS>0) - (Tr_s_tmp - Tr_s_hycom(t))/sum(sum(areaS));
obcs_corr(:,:,t) = tmp;

% calc transport from u_corr
Tr_n(t) = nansum(nansum(obcn_corr(:,:,t).*areaN));
Tr_s(t) = nansum(nansum(obcs_corr(:,:,t).*areaS));

end
end


if strcmp(param,'U')
for t=1:nt_obc

maskW=squeeze(hFacW(2,:,:));
areaW = (DYG(2,:)'*DRF').*maskW;
Tr_w_tmp = nansum(nansum(obcw_reg(:,:,t).*areaW));
tmp = obcw_reg(:,:,t);
tmp(maskW>0) = tmp(maskW>0) - (Tr_w_tmp - Tr_w_hycom(t))/sum(sum(areaW));
obcw_corr(:,:,t) = tmp;
Tr_w(t) = nansum(nansum(obcw_corr(:,:,t).*areaW));
Tr_tot(t) = Tr_s(t) + Tr_w(t) - Tr_n(t);
tmp = obcw_corr(:,:,t);
tmp(maskW>0) = tmp(maskW>0) - Tr_tot(t)/sum(sum(areaW));
obcw_zeronet(:,:,t) = tmp;

end

% match AVISO

% 2) load AVISO

load /data/SO6/TPOSE/constraints/aviso/aviso_tot_MSLA_eqpac.mat
time_aviso = tsave; clear tsave
lon_aviso = xsave; clear xsave
lat_aviso = ysave; clear ysave
ssh_aviso = permute(dsave,[3 2 1]); clear dsave
ssh_aviso(ssh_aviso>1e18)=NaN;
nt_aviso = length(time_aviso);

if 0
% not needed until 2020/6/1
% near real time 
load /data/averdy/tpose/constraints/SSH/AVISO/aviso_tot_MSLA_nrt_eqpac
ind = find(tsave==time_aviso(end)+1);
tsave = tsave(ind:end);
dsave = dsave(ind:end,:,:);
nt2 = length(tsave);
time_aviso = [time_aviso(:); tsave(:)]; clear tsave
ssh_aviso(:,:,nt_aviso+1:nt_aviso+nt2) = permute(dsave,[3 2 1]); clear dsave
end % if 0


% 3) select time period
% here: starting 1/1/2019

start_obcs = my_obc_date(1);
ind = find(time_aviso==start_obcs-3);
ssh_aviso = ssh_aviso(:,:,ind:5:end);
time_aviso = time_aviso(ind:5:end);
nt_aviso = length(time_aviso);
if nt_obc<nt_aviso
 time_aviso = time_aviso(1:nt_obc+1);
 ssh_aviso = ssh_aviso(:,:,1:nt_obc+1);
end
nt_aviso = length(time_aviso);


% 4) interpolate to model grid

ssh_aviso_interp = zeros(nx_obc,ny_obc,nt_aviso);
for t=1:nt_aviso
  ssh_aviso_interp(:,:,t) = interp2(lon_aviso',lat_aviso,ssh_aviso(:,:,t)',XC(:,1),YC(1,:),'nearest')';
end


% 5) domain-mean SSH

area = repmat(RAC.*hFacC(:,:,1),[1 1 nt_aviso]);
ssh_mean = squeeze(nansum(nansum(ssh_aviso_interp.*area,1),2))/sum(sum(area(:,:,1)));
%figure;plot(time_aviso,ssh_mean*100);
%set(gca,'xtick',datenum(2010:2017,1,1),'xticklabel',2010:2017);
%title('SSH domain average, from AVISO (cm)');


% 6) calculate transport
% diff(SSH [m] x area [m^2]) / 5 days [s] = [m^3/s]

tr = diff(squeeze(nansum(nansum(ssh_aviso_interp.*area,1),2)))/(86400*5);


% 7) infer correction to velocity 
% corr [m/2] = tr [m^3/s] / obc area [m^2]
% apply correction only on wet points (where mask>0)

maskW = squeeze(hFacW(2,:,:));
areaW = (DYG(2,:)'*DRF').*maskW;
corr = tr/sum(sum(areaW));

%for t=1:nt_obc
for t=1:nt_aviso-1
 tmp = obcw_zeronet(:,:,t);
 tmp(maskW>0) = tmp(maskW>0) + corr(t);
 obcw_aviso(:,:,t) = tmp;
end

end % if U



if np>2

fid=fopen([param '_obcn_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcn_reg,'single');
fclose(fid);

fid=fopen([param '_obcs_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcs_reg,'single');
fclose(fid);

fid=fopen([param '_obcw_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_reg,'single');
fclose(fid);


elseif np==1

fid=fopen([param '_obcn_bal_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcn_corr,'single');
fclose(fid);

fid=fopen([param '_obcs_bal_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcs_corr,'single');
fclose(fid);

fid=fopen([param '_obcw_zeros_ccs_' datestr '.bin'],'w','b');
fwrite(fid,zeros(ny_obc,nz_obc,nt_obc),'single');
fclose(fid);

elseif np==2

fid=fopen([param '_obcw_bal_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_corr,'single');
fclose(fid);

fid=fopen([param '_obcw_zeronet_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_zeronet,'single');
fclose(fid);

if 1
fid= fopen([param '_obcw_aviso_hycom_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_aviso,'single');
fclose(fid);
end

fid=fopen([param '_obcn_zeros_ccs_' datestr '.bin'],'w','b');
fwrite(fid,zeros(nx_obc,nz_obc,nt_obc),'single');
fclose(fid);

fid=fopen([param '_obcs_zeros_ccs_' datestr '.bin'],'w','b');
fwrite(fid,zeros(nx_obc,nz_obc,nt_obc),'single');
fclose(fid);

end


end % np

netcdf.close(ncid);



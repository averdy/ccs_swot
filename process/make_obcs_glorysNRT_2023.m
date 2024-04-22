% make OBCS by interpolating glorys


datestr = 'april2023';

addpath ~/scripts_m

% read model grid

load /data/SO6/CCS/grid/grid XC YC RC DXG DYG DRF hFacS hFacW hFacC RAC
[nx_obc ny_obc nz_obc]=size(hFacC);


cd('/data/SO6/CCS/input_processing/CMEMS_forecast/');
ncloads('cmems_t_obcn_april2023.nc');
clear thetao
nx=length(longitude);
ncloads('cmems_t_obcw_april2023.nc','latitude');
ny=length(latitude);
longitude=longitude+360;
depth=[-depth; -6000]; 

time = time/24+datenum(1950,1,1);

my_obc_time = datenum(2023,4,1):5:datenum(2023,4,30);
nt_obc = length(my_obc_time);


% loop parameters

params = {'v', 'u', 's', 'theta'};
for np = 1:length(params)

clear obc*

param = params{np};

if np>2
obcn_all = ncread(['cmems_' param(1) '_obcn_april2023.nc'],[param 'o']);
obcs_all = ncread(['cmems_' param(1) '_obcs_april2023.nc'],[param 'o']);
obcw_all = ncread(['cmems_' param(1) '_obcw_april2023.nc'],[param 'o']);
else
obcn_all = ncread(['cmems_uv_obcn_april2023.nc'],[param 'o']);
obcs_all = ncread(['cmems_uv_obcs_april2023.nc'],[param 'o']);
obcw_all = ncread(['cmems_uv_obcw_april2023.nc'],[param 'o']);
end

obcn_all=squeeze(obcn_all);
obcs_all=squeeze(obcs_all);
obcw_all=squeeze(obcw_all);

if np<4
obcn_all=obcn_all(1:nx,:,:);
obcs_all=obcs_all(1:nx,:,:);
end

nz=size(obcn_all,2);

% repeat bottom layer
obcn_all(:,nz+1,:)=obcn_all(:,nz,:);
obcs_all(:,nz+1,:)=obcs_all(:,nz,:);
obcw_all(:,nz+1,:)=obcw_all(:,nz,:);


cnt=1;
for t=1:nt_obc %-1
 % read 2 days before and after to make 5-day average
 % except first time
 if t==1
  ind=cnt:cnt+2;
 elseif t==nt_obc
  ind=cnt-2:cnt+1; 
 else
  ind=cnt-2:cnt+2;
 end
cnt=cnt+5;

obcn(:,:,t)=nanmean(obcn_all(:,:,ind),3);
obcs(:,:,t)=nanmean(obcs_all(:,:,ind),3);
obcw(:,:,t)=nanmean(obcw_all(:,:,ind),3);



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
tmp = interp2(longitude,depth,obcs_tmp(:,:,t)',XC(:,1)',RC)';
tmp(isnan(tmp))=0;			
obcs_reg(:,:,t) = tmp.*maskS;
tmp_s=tmp_s+sum(sum(obcs_reg(:,:,t)==0&maskS>0));			
tmp = interp2(longitude,depth,obcn_tmp(:,:,t)',XC(:,1)',RC)';
tmp(isnan(tmp))=0;			
obcn_reg(:,:,t) = tmp.*maskN;
tmp_n=tmp_n+sum(sum(obcn_reg(:,:,t)==0&maskN>0));
tmp = interp2(latitude,depth,obcw_tmp(:,:,t)',YC(1,:),RC)';
tmp(isnan(tmp))=0;			
obcw_reg(:,:,t) = tmp.*maskW;
tmp_w=tmp_w+sum(sum(obcw_reg(:,:,t)==0&maskW>0));
end
tmp_s,tmp_n,tmp_w


if strcmp(param,'v')
for t=1:nt_obc

maskN = squeeze(hFacS(:,end,:));
maskS = squeeze(hFacS(:,2,:));
areaN = ((DXG(:,end)*DRF').*maskN);
areaS = ((DXG(:,2)*DRF').*maskS);
Tr_n(t) = nansum(nansum(obcn_reg(:,:,t).*areaN));
Tr_s(t) = nansum(nansum(obcs_reg(:,:,t).*areaS));

end
end


if strcmp(param,'u')
for t=1:nt_obc

maskW=squeeze(hFacW(2,:,:));
areaW = (DYG(2,:)'*DRF').*maskW;
Tr_w(t) = nansum(nansum(obcw_reg(:,:,t).*areaW));
Tr_tot(t) = Tr_s(t) + Tr_w(t) - Tr_n(t);
tmp = obcw_reg(:,:,t);
tmp(maskW>0) = tmp(maskW>0) - Tr_tot(t)/sum(sum(areaW));
obcw_zeronet(:,:,t) = tmp;

end

if 0
figure;plot(Tr_n);hold on;plot(Tr_s);plot(Tr_w);plot(Tr_tot,'k');
for t=1:nt_obc
Tr_w(t) = nansum(nansum(obcw_zeronet(:,:,t).*areaW));
Tr_tot(t) = Tr_s(t) + Tr_w(t) - Tr_n(t);
end
figure;plot(Tr_n);hold on;plot(Tr_s);plot(Tr_w);plot(Tr_tot,'k');
end


if 0
% match AVISO

% 2) load AVISO

if 1
load /data/SO6/CCS/constraints/aviso/aviso_tot_MSLA_ccs.mat
time_aviso = tsave; clear tsave
lon_aviso = xsave; clear xsave
lat_aviso = ysave; clear ysave
ssh_aviso = permute(dsave,[3 2 1]); clear dsave
ssh_aviso(ssh_aviso>1e18)=NaN;
nt_aviso = length(time_aviso);
end

if 0
% near real time 
load /data/SO6/CCS/constraints/aviso/aviso_tot_MSLA_nrt_ccs
ind = find(tsave==time_aviso(end)+1);
tsave = tsave(ind:end);
dsave = dsave(ind:end,:,:);
nt2 = length(tsave);
time_aviso = [time_aviso(:); tsave(:)]; clear tsave
ssh_aviso(:,:,nt_aviso+1:nt_aviso+nt2) = permute(dsave,[3 2 1]); clear dsave
end % if 0

if 0
load /data/SO6/CCS/constraints/aviso/aviso_tot_MSLA_nrt_ccs.mat
time_aviso = tsave; clear tsave
lon_aviso = xsave; clear xsave
lat_aviso = ysave; clear ysave
ssh_aviso = permute(dsave,[3 2 1]); clear dsave
ssh_aviso(ssh_aviso>1e18)=NaN;
nt_aviso = length(time_aviso);
end % if 0


% 3) select time period
% here: starting 1/1/2016

start_obcs = my_obc_date(1);
%ind = find(time_aviso==start_obcs-3);
ind = find(time_aviso==start_obcs);
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
%tr_tot = tr + Tr_s - Tr_n;

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

end % if 0

end % if U




if 0
tmp = obcn_reg;
newobc = tmp(:,:,1);
newobc(:,:,2:nt_obc+1) = tmp;
obcn_reg = newobc;

tmp = obcs_reg;
newobc = tmp(:,:,1);
newobc(:,:,2:nt_obc+1) = tmp;
obcs_reg = newobc;

tmp = obcw_reg;
newobc = tmp(:,:,1);
newobc(:,:,2:nt_obc+1) = tmp;
obcw_reg = newobc;
end % if 0



fid=fopen([param(1) '_obcn_glorysNRT_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcn_reg,'single');
fclose(fid);

fid=fopen([param(1) '_obcs_glorysNRT_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcs_reg,'single');
fclose(fid);

fid=fopen([param(1) '_obcw_glorysNRT_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_reg,'single');
fclose(fid);


if np==2

fid=fopen([param(1) '_obcw_zeronet_glorysNRT_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_zeronet,'single');
fclose(fid);

if 0
fid= fopen([param(1) '_obcw_aviso_glorysNRT_ccs_' datestr '.bin'],'w','b');
fwrite(fid,obcw_aviso,'single');
fclose(fid);
end


end


end % np




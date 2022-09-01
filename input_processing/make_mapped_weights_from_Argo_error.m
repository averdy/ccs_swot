clear all 
close all

addpath ~/scripts_m/

% read maps 2004-2018
ncloads('/data/averdy/datasets/ARGO/RG_ArgoClim_Temperature_2019.nc');
ncloads('/data/averdy/datasets/ARGO/RG_ArgoClim_Salinity_2019.nc');
% convert pressure to depth
for y=1:130
 z(:,y,:)=repmat(sw_dpth(PRESSURE,LATITUDE(y)),[1 1 360]);
end

yrmin=2010; yrmax=2018;
nt=(yrmax-yrmin+1)*12;

% grid on which we bin for error calculation
xbox=LONGITUDE; ybox=LATITUDE;
% xbox=25:10:380; ybox=-64:2:64;

% initialize gridded fields
nx=length(xbox); ny=length(ybox); 
nz=97; % number of levels for binning for error calculation
% number of profiles per grid box
nprofT = zeros(nx,ny,nz);
nprofS = zeros(nx,ny,nz);
sum_misfitT_sq = zeros(nx,ny,nz);
sum_misfitS_sq = zeros(nx,ny,nz);

tcount=0;

for yr=2010:2018
 display(yr)
 % read profiles
 ncloads(['/data/shared_data/Interp_Profiles/world/assim/PFL/USGO_CA_' num2str(yr) '_PFL_D.nc']);
 prof_T(prof_Tweight==0)=NaN;
 prof_S(prof_Sweight==0)=NaN;
 prof_T(prof_T<0)=NaN;
 prof_S(prof_S<0)=NaN;
 for mon=1:12
 tcount=tcount+1;
 display(tcount)
 % select profs
 ind=find(prof_date>datenum(yr,mon,0) & prof_date<datenum(yr,mon+1,0));
 prof_Tm=prof_T(ind,:);
 prof_Sm=prof_S(ind,:);
 prof_latm=prof_lat(ind);
 prof_lonm=prof_lon(ind);
  t = (yr-2004)*12+mon;
  T = ARGO_TEMPERATURE_MEAN+squeeze(ARGO_TEMPERATURE_ANOMALY(t,:,:,:));
  S = ARGO_SALINITY_MEAN+squeeze(ARGO_SALINITY_ANOMALY(t,:,:,:));
  T(T<0)=NaN;
  S(S<0)=NaN;
  % convert to potential temperature
  PR=0;
  for x=1:360
   THETA(:,:,x) = sw_ptmp(S(:,:,x),T(:,:,x),PRESSURE,PR);
  end
  % calculate misfits
  np=length(prof_latm);
  misfitT_cur_month=NaN(nx,ny,nz,np);
  misfitS_cur_month=NaN(nx,ny,nz,np);
  % subsample maps
  clear T_map_equi S_map_equi
  for n=1:np
   [~,i] = min(abs(prof_lonm(n)-LONGITUDE));
   [~,j] = min(abs(prof_latm(n)-LATITUDE));
   T_map_equi = interp1(z(:,j,i),THETA(:,j,i),prof_depth);
   S_map_equi = interp1(z(:,j,i),S(:,j,i),prof_depth);
   misfitT = T_map_equi-prof_Tm(n,:)';
   misfitS = S_map_equi-prof_Sm(n,:)';
   [~,y]=min(abs(prof_latm(n)-ybox));
   [~,x]=min(abs(prof_lonm(n)-xbox));
   if sum(~isnan(misfitT))>0
    nprofT(x,y,:) = squeeze(nprofT(x,y,:))+(~isnan(misfitT));
    misfitT(isnan(misfitT)) = 0;
    sum_misfitT_sq(x,y,:) = squeeze(sum_misfitT_sq(x,y,:))+misfitT.^2;
   end
   if sum(~isnan(misfitS))>0
    nprofS(x,y,:) = squeeze(nprofS(x,y,:))+(~isnan(misfitS));
    misfitS(isnan(misfitS)) = 0;
    sum_misfitS_sq(x,y,:) = squeeze(sum_misfitS_sq(x,y,:))+misfitS.^2;
   end
  end % profs
 end % month
end % year
% rms misfits
rms_misfitT = sqrt(sum_misfitT_sq./nprofT);
rms_misfitS = sqrt(sum_misfitS_sq./nprofS);

%save('/data/averdy/tpose/weights/rms_error_2018.mat','sum_misfitT_sq','sum_misfitS_sq','xbox','ybox');
save('/data/SO6/CCS/Argo_rms_error.mat','rms_misfitT','rms_misfitS','xbox','ybox');


% keep only well populated bins
nthresh = 50;
rms_misfitT(nprofT<nthresh) = NaN;
rms_misfitS(nprofS<nthresh) = NaN;

% vertical profile
rms_misfitT_vs_z = squeeze(sqrt(nansum(nansum(sum_misfitT_sq,2),1)./sum(sum(nprofT,2),1)));
rms_misfitS_vs_z = squeeze(sqrt(nansum(nansum(sum_misfitS_sq,2),1)./sum(sum(nprofS,2),1)));

% optimal interpolation
load /data/SO6/CCS/grid/grid XC YC RC Depth
[xx yy]=meshgrid(xbox,ybox);
xx=reshape(xx,[1,nx*ny]);
yy=reshape(yy,[1,nx*ny]);
lat = 0;
sigman = 0.3;
icovar = 1;
dcors = 1500;

for k=2:nz
 icovar = rms_misfitT_vs_z(k);
 sigman = icovar/2;
 tmp=reshape(rms_misfitT(:,:,k)',[1,nx*ny]);
 [that2tmp Phat2] = CheapObjMap(xx,yy,tmp,lat,sigman,icovar,dcors);
 that2(:,:,k) = reshape(that2tmp,ny,nx);
 that2(that2<0.05) = 0.05;
 tmp=reshape(rms_misfitS(:,:,k)',[1,nx*ny]);
 [shat2tmp Phat2] = CheapObjMap(xx,yy,tmp,lat,sigman,icovar,dcors);
 shat2(:,:,k) = reshape(shat2tmp,ny,nx);
 shat2(shat2<0.05) = 0.05;
 % interp tpose domain
 error_Ti=interp2(xbox,ybox,that2(:,:,k),XC(:,1),YC(1,:))';
 error_Ti(Depth==0)=0;
 error_T(:,:,k) = error_Ti;
 error_Si=interp2(xbox,ybox,shat2(:,:,k),XC(:,1),YC(1,:))';
 error_Si(Depth==0)=0;
 error_S(:,:,k) = error_Si;
end
error_T(:,:,1) = error_T(:,:,2);
error_S(:,:,1) = error_S(:,:,2);






ST
fid=fopen([opath 'SST_fromArgo.err'],'w','b');;;;;
fwrite(fid,error_T(:,:,1),'single');;
fclose(fid);

% interpolate vertically
error_T = reshape(error_T,nx*ny,25);
error_T_interp = interp1(-depth,error_T',RC)';
error_T_interp = reshape(error_T_interp,564,168,51);

% fill in upper levels
error_T_interp(:,:,1) = error_T_interp(:,:,3);
error_T_interp(:,:,2) = error_T_interp(:,:,3);

% multiply by 10
error_T = error_T_interp*10;

% save
fid=fopen([opath 'ic_T_fromArgo_x10.err'],'w','b');
fwrite(fid,error_T,'single');
fclose(fid);

load /data/averdy/tpose/weights/error_S_fromArgo.mat

% interpolate vertically
error_S = reshape(error_S,nx*ny,25);
error_S_interp = interp1(-depth,error_S',RC)';
error_S_interp = reshape(error_S_interp,564,168,51);

% fill in upper levels
error_S_interp(:,:,1) = error_S_interp(:,:,3);
error_S_interp(:,:,2) = error_S_interp(:,:,3);

% multiply by 10
error_S = error_S_interp*10;

% save
fid=fopen([opath 'ic_S_fromArgo_x10.err'],'w','b');
fwrite(fid,error_S,'single');
fclose(fid);


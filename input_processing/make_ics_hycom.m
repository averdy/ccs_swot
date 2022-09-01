% read HYCOM 
% snapshop for an arbitrary day
% interpolate to model grid


dates_to_read = datenum(2019,10,1);
nt = length(dates_to_read);

for d=1:nt

 newhycom = 0;

 if dates_to_read(d)<=datenum(2011,2,1)
  expt = 'expt_90.8';
 elseif dates_to_read(d)>datenum(2011,2,1) & dates_to_read(d)<=datenum(2013,8,20)
  expt = 'expt_90.9';
 elseif dates_to_read(d)>datenum(2013,8,20) & dates_to_read(d)<=datenum(2014,4,4)
  expt = 'expt_91.0';
 elseif dates_to_read(d)>datenum(2014,4,4) & dates_to_read(d)<=datenum(2016,4,17)
  expt = 'expt_91.1';
 elseif dates_to_read(d)>datenum(2016,4,17) & dates_to_read(d)<=datenum(2018,11,27)
  expt = 'expt_91.2';
 elseif dates_to_read(d)>datenum(2018,11,27) 
  expt = 'expt_93.0';
  newhycom = 1;
 end

% open data file
if newhycom==0
 OpenDAP_URL = ['http://tds.hycom.org/thredds/dodsC/GLBa0.08/' expt];
else
 OpenDAP_URL = ['http://tds.hycom.org/thredds/dodsC/GLBv0.08/' expt];
end
ncid = netcdf.open(OpenDAP_URL,'NOWRITE');
prec = 'double';

% Variable id in data set
if newhycom==0
nc_depth = 2;
nc_date = 4;
nc_lat = 5;
nc_lon = 6;
nc_ssh = 11;
nc_salinity = 14; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,14)
nc_temperature = 15; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,15)
nc_u = 16; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,16)
nc_v = 17; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,17)
else
nc_depth = 4;
nc_date = 2;
nc_lat = 0;
nc_lon = 1;
nc_ssh = 5;
nc_salinity = 8; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,8)
nc_temperature = 6; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,6)
nc_u = 10; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,10)
nc_v = 12; % [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,12)
end


if newhycom==0
else
 % 230E to 246.2E;
 xmin = 2875-1; xmax = 3079-1; nx = xmax-xmin+1;
 % 27N to 43.1N
 ymin = 1838-1; ymax = 2079-1; ny = ymax-ymin+1;
 zmin = 0; zmax = 39; nz = zmax-zmin+1;
end

hycom_depth = netcdf.getVar(ncid,nc_depth,prec);
hycom_depth = -1*hycom_depth;
% HYCOM GOES ONLY UPTO 5500 m, so adding one more layer with same field as
hycom_depth(end+1) = -6500;

% middle of layer: 
hycom_depth = (hycom_depth(1:end-1) + hycom_depth(2:end))/2;

hycom_lat = netcdf.getVar(ncid,nc_lat,ymin,ny,prec);

hycom_lon = netcdf.getVar(ncid,nc_lon,xmin,nx,prec);
if hycom_lon < 0; hycom_lon = hycom_lon + 360; end


hycom_date = netcdf.getVar(ncid,nc_date);
if newhycom==0
 hycom_date = datenum(datevec(num2str(hycom_date),'yyyymmdd'));
else
 % units are hours since 2000-1-1
 hycom_date = datenum(2000,1,1)+hycom_date/24;
end

if 0
if newhycom==0
display(['first day: ' num2str(hycom_date(1))])
display(['last day: ' num2str(hycom_date(end))])
else
display(['first day: ' datevec(datenum(2000,0,0)+hycom_date(1)/24)])
display(['last day: ' datevec(datenum(2000,0,0)+hycom_date(end)/24)])
end
display(['nt=' num2str(nt) ' , nthycom=' num2str(nthycom)])
end

ind = find(dates_to_read(d)==hycom_date);
if ind>0
 T = (netcdf.getVar(ncid,nc_temperature,[xmin,ymin,zmin,ind-1],[nx,ny,nz,1],prec)); 
 S = (netcdf.getVar(ncid,nc_salinity,[xmin,ymin,zmin,ind-1],[nx,ny,nz,1],prec)); 
 U = (netcdf.getVar(ncid,nc_u,[xmin,ymin,zmin,ind-1],[nx,ny,nz,1],prec)); 
 V = (netcdf.getVar(ncid,nc_v,[xmin,ymin,zmin,ind-1],[nx,ny,nz,1],prec)); 
else
 T = NaN(nx,ny,nz);
 S = NaN(nx,ny,nz);
 U = NaN(nx,ny,nz);
 V = NaN(nx,ny,nz);
end


% scaling
% from Ganesh
scale = [0.001, 0.001, 0.001, 0.001];
offset = [20, 20, 0, 0];


miss_value = -3e4;
T(T <= miss_value) = nan;
S(S <= miss_value) = nan;
U(U <= miss_value) = nan;
V(V <= miss_value) = nan;

T=T*scale(1)+offset(1);
S=S*scale(2)+offset(2);
U=U*scale(3)+offset(3);
V=V*scale(4)+offset(4);


load /data/SO6/CCS/grid/grid.mat XC YC RC hFacC
[nx1,ny1] = size(XC); nz1=length(RC);


% first interpolate horizontally
for k=1:nz
 Ttmp = T(:,:,k);
 Stmp = S(:,:,k);
 Utmp = U(:,:,k);
 Vtmp = V(:,:,k);

         % prepare for interpolation
            for j=ny-1:-1:1
            for i=2:nx-1
        if(isnan(Ttmp(i,j))) 
              Ttmp(i,j) = nanmean(Ttmp(i-1:i+1,j));
              Stmp(i,j) = nanmean(Stmp(i-1:i+1,j));
              Utmp(i,j) = nanmean(Utmp(i-1:i+1,j));
              Vtmp(i,j) = nanmean(Vtmp(i-1:i+1,j));
        end
        if(isnan(Ttmp(i,j))) 
              Ttmp(i,j) = Ttmp(i,j+1);
              Stmp(i,j) = Stmp(i,j+1);
              Utmp(i,j) = Utmp(i,j+1);
              Vtmp(i,j) = Vtmp(i,j+1);
        end
            end
            end

 Ti(:,:,k) = interp2(hycom_lat,hycom_lon,Ttmp,YC,XC);
 Si(:,:,k) = interp2(hycom_lat,hycom_lon,Stmp,YC,XC);
 Ui(:,:,k) = interp2(hycom_lat,hycom_lon,Utmp,YC,XC);
 Vi(:,:,k) = interp2(hycom_lat,hycom_lon,Vtmp,YC,XC);
end

    % Then interpolate vertically
        for k=2:nz
            for j=1:ny1
            for i=1:nx1
            if(isnan(Ti(i,j,k))) 
            Ti(i,j,k) = Ti(i,j,k-1);
            Si(i,j,k) = Si(i,j,k-1);
            Ui(i,j,k) = Ui(i,j,k-1);
            Vi(i,j,k) = Vi(i,j,k-1);
        end
            end
            end    
        end
            
    Ti = reshape(Ti,nx1*ny1,nz)';
    Ti = interp1(hycom_depth,Ti,RC,'linear','extrap');
    Ti = reshape(Ti',nx1,ny1,nz1);

    Si = reshape(Si,nx1*ny1,nz)';
    Si = interp1(hycom_depth,Si,RC,'linear','extrap');
    Si = reshape(Si',nx1,ny1,nz1);

    Ui = reshape(Ui,nx1*ny1,nz)';
    Ui = interp1(hycom_depth,Ui,RC,'linear','extrap');
    Ui = reshape(Ui',nx1,ny1,nz1);

    Vi = reshape(Vi,nx1*ny1,nz)';
    Vi = interp1(hycom_depth,Vi,RC,'linear','extrap');
    Vi = reshape(Vi',nx1,ny1,nz1);


Ti(hFacC==0)=0;
Si(hFacC==0)=0;
Ui(hFacC==0)=0;
Vi(hFacC==0)=0;


fid=fopen(['/data/SO6/CCS/ics/T_ccs_1oct2019.bin'],'w','b');
fwrite(fid,Ti,'single');
fclose(fid);

fid=fopen(['/data/SO6/CCS/ics/S_ccs_1oct2019.bin'],'w','b');
fwrite(fid,Si,'single');
fclose(fid);

fid=fopen(['/data/SO6/CCS/ics/U_ccs_1oct2019.bin'],'w','b');
fwrite(fid,Ui,'single');
fclose(fid);

fid=fopen(['/data/SO6/CCS/ics/V_ccs_1oct2019.bin'],'w','b');
fwrite(fid,Vi,'single');
fclose(fid);



end % for d

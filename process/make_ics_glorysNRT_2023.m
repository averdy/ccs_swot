% read GLORYS
% interpolate to model grid



cd /data/SO6/CCS/input_processing/glorys/CMEMS_forecast_data/
addpath ~/scripts_m
ncloads('cmems_t_10may2023.nc');
thetao=ncread('cmems_t_10may2023.nc','thetao');
so=ncread('cmems_s_10may2023.nc','so');
uo=ncread('cmems_uv_10may2023.nc','uo');
vo=ncread('cmems_uv_10may2023.nc','vo');
[nx,ny,nz]=size(thetao);

depth(end)=6000;
longitude=longitude+360;


load /data/SO6/CCS/grid/grid.mat XC YC RC hFacC
[nx1,ny1] = size(XC); nz1=length(RC);


% first interpolate horizontally
for k=1:nz-1
 Ttmp = thetao(:,:,k);
 Stmp = so(:,:,k);
 Utmp = uo(:,:,k);
 Vtmp = vo(:,:,k);

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

 Ti(:,:,k) = interp2(latitude,longitude,Ttmp,YC,XC);
 Si(:,:,k) = interp2(latitude,longitude,Stmp,YC,XC);
 Ui(:,:,k) = interp2(latitude,longitude,Utmp,YC,XC);
 Vi(:,:,k) = interp2(latitude,longitude,Vtmp,YC,XC);
end

Ti(:,:,nz)=Ti(:,:,nz-1);
Si(:,:,nz)=Si(:,:,nz-1);
Ui(:,:,nz)=Ui(:,:,nz-1);
Vi(:,:,nz)=Vi(:,:,nz-1);


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
    Ti = interp1(-depth,Ti,RC,'linear','extrap');
    Ti = reshape(Ti',nx1,ny1,nz1);

    Si = reshape(Si,nx1*ny1,nz)';
    Si = interp1(-depth,Si,RC,'linear','extrap');
    Si = reshape(Si',nx1,ny1,nz1);

    Ui = reshape(Ui,nx1*ny1,nz)';
    Ui = interp1(-depth,Ui,RC,'linear','extrap');
    Ui = reshape(Ui',nx1,ny1,nz1);

    Vi = reshape(Vi,nx1*ny1,nz)';
    Vi = interp1(-depth,Vi,RC,'linear','extrap');
    Vi = reshape(Vi',nx1,ny1,nz1);


Ti(hFacC==0)=0;
Si(hFacC==0)=0;
Ui(hFacC==0)=0;
Vi(hFacC==0)=0;


fid=fopen(['/data/SO6/CCS/linked_files/may2023/T_ccs_glorysNRT_10may2023.bin'],'w','b');
fwrite(fid,Ti,'single');
fclose(fid);

fid=fopen(['/data/SO6/CCS/linked_files/may2023/S_ccs_glorysNRT_10may2023.bin'],'w','b');
fwrite(fid,Si,'single');
fclose(fid);

fid=fopen(['/data/SO6/CCS/linked_files/may2023/U_ccs_glorysNRT_10may2023.bin'],'w','b');
fwrite(fid,Ui,'single');
fclose(fid);

fid=fopen(['/data/SO6/CCS/linked_files/may2023/V_ccs_glorysNRT_10may2023.bin'],'w','b');
fwrite(fid,Vi,'single');
fclose(fid);




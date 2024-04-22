clear all

%WHICH GRID
x = ncread('/data/soccom/SO3input/RUNOFF/DATA/runoff.15JUNE2009.nc','lon');
y = ncread('/data/soccom/SO3input/RUNOFF/DATA/runoff.15JUNE2009.nc','lat');
F = ncread('/data/soccom/SO3input/RUNOFF/DATA/runoff.15JUNE2009.nc','Foxx_o_roff');%(kg/s)/m^2
F = F./1000; %m/s
x(x<0)=x(x<0)+360;
F(F>9999)=0;
F(isnan(F))=0;
%MAKE CLIM OUT OF LAST 30 YEARS
%MAP
  fld_out = 'runoff_core_cnyf2p0_ccs_ms.bin';
  load /data/SO6/CCS/grid/grid hFacC XC YC
  OBC = 6;
[NX NY NZ] = size(hFacC);
xi = XC(:,1);yi = YC(1,:);
FLOW=zeros(NX,NY);
%IS SMOOTH FILE SO WE WILL JUST USE INTERP
FLOW = interp2(y,x,F,yi,xi,'linear');
FLOW(1:OBC,:) = 0;
FLOW(:,NY-OBC+1:NY) = 0;
FLOW(:,1:OBC) = 0;

FLOW(isnan(FLOW))=0;
max(FLOW(:))  %want less than 1e-6
sum(FLOW(:))

figure(2);clf; 
pcolor(xi,yi,FLOW'); shading flat
caxis([0 1e-7]); colorbar
 
% make 12 months
FLOW=repmat(FLOW,[1 1 12]);
 
fid=fopen([fld_out],'w','b');
fwrite(fid,FLOW,'single');
fclose(fid);


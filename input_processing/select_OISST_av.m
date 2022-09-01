% select_OISST.m
%
% function [Dsst,Derror,Dmask,DLon,Dlat] =  select_OISST(filename,domain)
%   EXAMPLE:  [Dsst,Derror,Dmask,DLon,Dlat] =  select_OISST('data/tmi_amsre.fusion.2002.152.v03.gz','SO');
%
% matlab program (written for R2009b) to call the OISST provided matlab program
%     read_rss_oisst.m
%
% and to return sst, error, mask for the specified SIO domain
%
%  input:  oisst_fn = the oist file name to read, may include a path. example: tmi_amsre.fusion.2002.152.v03.gz
%  input:  domain = the domain name of interest.  example: 'SO' for southern ocean, or 'none'
%		to use the area specified by editing this file.
%		Currently,  only 'SO' and 'none' are accepted, but 'none' exits the program
%
%  returns sst, error, mask arrays for the domain of interest
% 	The MW only data files consist of 3 arrays sized 1440x720 (0.25 deg grid)
%
%  the grid description is NOT in the data file.  
%
%  Below from 
%     webpage http://www.remss.com/measurements/sea-surface-temperature/oisst-description
%    "center of the first cell of the 1440 column and 720 row map is at 0.125 E longitude and -89.875 latitude. 
%     The center of the second cell is 0.375 E longitude, -89.875 latitude."
%
%  The above grid description indicates that this is NOT a 25km grid, as the website states,  but a 0.25 degree grid.
%  See read_rss_oisst.m for further details of the data file.
%  
%
% The read_rss_oisst function returns:
% sst = array of real SST values
% error = array of interpolation error 
% mask = array data mask
%
%
% The mask should be used to identify sea ice, land, bad data.  
% Mask values are:
% Bit 1 = land
% Bit 2 = ice
% Bit 3 = IR data used
% Bit 4 = MW data used
% Bit 5 = bad data (near land)
%
% In this sample program the mask values are used to set data in the SST and Error arrays to:
%      252 = sea ice
%      254 = missing data, near land
%      255 = land mass
%


function [Dsst,Derror,Dmask,DLon,DLat] =  select_OISST_av(filename,domain)

if ~exist(filename)
   error = ['file not found: ', filename]
   sst=[];error=[];mask=[];
   return;
end;

% addpath download_example_programs/matlab

% read the filename
unz_filename = strrep(filename,'.gz','');
eval(['! gunzip ' filename ]);
[sst,error,mask] = read_rss_oisst(unz_filename);
eval(['! gzip ' unz_filename ]);


% To create a TEST OUTPUT IMAGE use:
if 0
str={'Land','Ice','IR data','MW data','Bad data'};
figure(1); subplot(2,1,1),imagesc(sst',[-3 34]); subplot(2,1,2),imagesc(sst',[250 255]);
figure(2);for i=1:5,ilnd=bitget(mask,i);subplot(1,5,i),imagesc(ilnd',[0 1]);title(str{i});end;
figure(3);subplot(2,1,1),imagesc(error',[0 1.1]);subplot(2,1,2),imagesc(error',[250 255]);
figure(4);imagesc(sst',[-3 34]); set(gca,'YDir','normal');
end


% To select the DOMAIN
% 
if  strmatch(lower(domain), 'none', 'exact')
%        user defined supplied here for  ymin,ymax,xmin,xmax
   disp('NOT programmed.  Edit  select_OISST.m and retry.  Quitting')
elseif  strmatch(domain, 'SO', 'exact')
   Dminy = -90;
   Dmaxy = -24.8;
   Dminx = 0;
   Dmaxx = 360; 
elseif  strmatch(domain, 'troppac', 'exact')
   Dminy = -27;
   Dmaxy = 31;
   Dminx = 103;
   Dmaxx = 293;
elseif  strmatch(domain, 'ccs', 'exact')
   Dminy = 26;
   Dmaxy = 44;
   Dminx = 229;
   Dmaxx = 247; 
else
   disp('Unknown domain. Only SO, ccs and troppac currently work')
end

% this from the comments in read_rss_oisst.m,  the oisst supplied program
%    "center of the first cell of the 1440 column and 720 row map is at 0.125 E longitude and -89.875 latitude. 
%     The center of the second cell is 0.375 E longitude, -89.875 latitude."
xmax = size(sst,1);
ymax = size(sst,2);
%xcell = 0.125:0.250:360; % grid cell values between 1 and xmax (where xmax = 4096 or 1440)
%ycell = -89.875:0.250:90; % grid cell values between 1 and ymax (where ymax = 2048 or 720)
xcell = 1:1440; % grid cell values between 1 and xmax (where xmax = 4096 or 1440)
ycell = 1:720; % grid cell values between 1 and ymax (where ymax = 2048 or 720)

dx=360./xmax;
dy=180./ymax;
CLon  = dx*xcell-dx/2.;% (degrees east), Center of grid cell Longitude
CLat  = dy*ycell-(90+dy/2.); % (-90 to 90); Center of grid cell Latitude

% NOTE min,max sst in OISST are -3,34
% figure(5); imagesc(CLon,CLat,sst',[-3 34]); set(gca,'YDir','normal');

% reset to return only domain areas:
iDCGX = find( (CLon>=Dminx) & (CLon<=Dmaxx)); % index of Domain_Center_grid_X
iDCGY = find( (CLat>=Dminy) & (CLat<=Dmaxy)); % index of Domain_Center_grid_Y

DLon = CLon(iDCGX); % Domain_Center_grid_X 
DLat = CLat(iDCGY); % Domain_Center_grid_Y

% figure(6); imagesc(DLon,DLat,sst(iDCGX,iDCGY)',[-3 34]); set(gca,'YDir','normal');


% return the Domain arrays
Dsst = sst(iDCGX,iDCGY);
Derror = error(iDCGX,iDCGY);
Dmask = mask(iDCGX,iDCGY);
% DLon, DLat set above

% figure(7); imagesc(DLon,DLat,Dsst',[-3 34]); set(gca,'YDir','normal');

return;



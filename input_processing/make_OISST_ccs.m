clear all 
close all

% path to data
dpath='/project_shared/OISST/';

addpath ~/scripts_m
load /data/SO6/CCS/grid/grid XC YC hFacC

for year= 2019 % 2010:2018

 if year==2016 | year==2012 | year==2008 | year==2004
  ndays = 366;
 else
  ndays = 365;
 end

 fid = fopen(['/data/SO6/CCS/constraints/SST/mw_fusion_ccs_' num2str(year)],'w','b');

 for day=1:ndays
  day_str = sprintf('%0.3i',day);
%  if ((year>2013 & day > 59) | (year>2014)) 
%   kind = 'rt';
%  else
   kind = 'v05.0';
%  end
  filename = [dpath 'data_v5/' num2str(year) '/mw.fusion.' num2str(year) '.' day_str '.' kind '.gz'];
  [Dsst,Derror,Dmask,DLon,Dlat] =  select_OISST_av(filename,'ccs');
  % missing values
  Dsst(Dsst>200) = NaN;
  SST = interp2(DLon,Dlat,Dsst',XC(:,1),YC(1,:))';
  SST(hFacC(:,:,1)==0)=-9999;
  SST(isnan(SST))=-9999;
  fwrite(fid,SST,'single');
 end

 fclose(fid);

end


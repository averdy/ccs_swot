% make weights for atm controls
% based on ERA5 variance

mypath = '/project_shared/ERA5/';

fields = {'dlw','dsw','rain','spfh2m','tmp2m_degC','u10m','v10m'};

wtmin = [1e-4 1e-5 1e14 1e5 0.05 0.01 0.01];
wtmax = [2e-2 4e-5 1e15 2e6 2 1.2 2];

load /data/SO6/CCS/grid/grid XC YC

mypath = '/project_shared/ERA5/';
outpath = '/data/SO6/CCS/';

% for smooth function:
addpath ~/scripts_m


lat(1) = -89.7848850;
dy = [0.2786766,0.2803199,0.2806917,0.2808328,0.2809011,0.2809392,0.2809626,0.2809781,0.2809888,0.2809965,0.2810023,0.2810067,0.2810102,0.2810129,0.2810151,ones(1,33)*0.2810258,ones(1,543)*0.2810302,ones(1,33)*0.2810258,0.2810151,0.2810129,0.2810102,0.2810067,0.2810023,0.2809965,0.2809888,0.2809781,0.2809626,0.2809392,0.2809011,0.2808328,0.2806917,0.2803199,0.2786766];
NY = length(dy)+1;
for i = 2:NY
  lat(i) = lat(i-1) + dy(i-1);
end
NX = 1280;
lon(1) = 0;
for i = 2:NX
  lon(i) = lon(i-1) + 0.2812500;
end

nx = length(lon); ny = length(lat);

% calc RMS of ERA atm state

for n=2:length(fields)

display(char(fields(n)))

f = char(fields(n));

cnt = 1;
%for yr=2010:2019
for yr=2010

display(num2str(yr))

fid=fopen([mypath 'ERA5_' f '_' num2str(yr)],'r','b');

if yr==2020 | yr==2016 | yr==2012 | yr==2008 | yr==2004
 nt=366*24;
else
 nt=365*24;
end

for t=1:nt

A = fread(fid,nx*ny,'single');
A = single(A);
A = reshape(A,nx,ny);

B(:,:,cnt) = interp2(lon,lat,A',XC(:,1),YC(1,:))';
cnt = cnt+1;

end

end % for yr

 Q = std(B,0,3);
 
 Q(Q==0)=min(Q(:));
 Q(isnan(Q))==min(Q(:));
 % S = Smooth2Dfnc(nx,ny,Q,1,[17 17],5,0);
 S=Q;

 wt = S.^-2;
 wt(wt<wtmin(n)) = wtmin(n);
 wt(wt>wtmax(n)) = wtmax(n);

 % figure; pcolor(wt'); shading flat; colorbar

fout = fopen([outpath 'ERA5_' char(fields(n)) '_wt_ccs.bin'],'w','b');
fwrite(fout,wt,'single');
fclose(fout);

fout = fopen([outpath 'ERA5_' char(fields(n)) '_std_ccs.bin'],'w','b');
fwrite(fout,Q,'single');
fclose(fout);

 clear Q S B wt 

end % for field





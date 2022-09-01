close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AVISO data are beeing updated by Caroline at /net/shared-data/
%% from 2013 May onwards. NO longer local copies
%% USING BRUC'es SCRIPT TO GET MAT FILE FOR GOM OR ANY PLACE
%% AVISO CHANGED DATA FORMAT IN APR 15, 2014,
%% THE NEW DATA FOR REGIONAL DOMAINS ARE SAVED IN AVISO_DATA_NEW/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opath = '/data/averdy/tpose/constraints/SSH/AVISO/';

%% start diary
fileOut = ['aviso_DT_region_pack_' datestr(date,29) '.log' ];
diary( [fileOut]);
disp('==========================');
disp(fileOut)

regions = {'eqpac'};nreg = length(regions);
for ar = 1:nreg
    reg_lab = regions{ar};
    disp(['region is ' reg_lab])
    
%    for rn = 1:2
    for rn=2
        if rn == 1
            madt = 1; msla = 0;
        elseif rn == 2
            madt = 0; msla = 1;
        end
        %% default database location
        dpath = ['/project_shared/DATA_DOWNLOADS/DATA/AVISO_DATA/AVISO_DATA_PHY/my.cmems-du.eu/Core/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/dataset-duacs-rep-global-merged-allsat-phy-l4/'];
        
        %% CRON JOB is crashed, so downloaded locally manually on nov 27, 2020
        %% MANUALLY DOWNLOADED on /project/ganeshgopal_hpc_data/AVISO_GG_20201126/NRT on NOV-27, 2020
        %% dpath = ['/project/ganeshgopal_hpc_data/AVISO_GG_20201126/DT/my.cmems-du.eu/Core/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/dataset-duacs-rep-global-merged-allsat-phy-l4/'];
        
        %% select for a region
        %% pick a region label
        %% will also be file label
        
        %% northwest pacific near taiwan
        %% reg_lab = 'nwpac';
        
        %% tropical pacific (+/- 10) near where p09 crosses the eq
        %% reg_lab = 'p09_eq';
        
        %% california current
        %% reg_lab = 'ccs';
        
        %% gulf of mexico
        %% reg_lab = 'gom';
        
        %% whole tropical pacific (huge if no subsample!)
        %% reg_lab = 'eqpac';
        
        %% okmc region
        %% reg_lab = 'okmc';
        
        %% BAY OF BENGAL
        %% reg_lab = 'bob'
        
        %% RED SEA
        %% reg_lab = 'red'
        
        %% INDIAN OCEAN
        %% reg_lab = 'ioc'
        
        %% HAWAII REGION
        %% reg_lab = 'hawaii'
        
        %% southern ocean 31S for Lynn
        %% reg_lab = 'SO_31S';
        
        %% whole tropical pacific (huge if no subsample!)
        %% reg_lab = 'bigeqpac';
        
        %% ceckr
        %% reg_lab = 'ceckr';
        
        %% default: don't subsample (x,y)
        i_sub_xy = [1 1];
        
        %% set limits based on choice of region
        if (strcmp(reg_lab,'eqpac'))
            disp('whole tropical pacific')
            x_min = 104;
            x_max = 292;
            y_min = -26.0;
            y_max =  30.0;
            %% sub-sample to keep file small (full domain is 500 Mb as of 5/26/10)
            %% at 1/4 degree, point count is: 60*4*150*4*52*20 = 149760000
            %% subsample to 1 degree
            %% i_sub_xy = [4 4];
            
        elseif (strcmp(reg_lab,'p09_eq'))
            disp('tropical pacific (+/- 10) near where p09 crosses the eq')
            x_min = 180;
            x_max = 220;
            y_min = -10.0;
            y_max =  10.0;
            
        elseif (strcmp(reg_lab,'okmc'))
            disp('okmc region: offshore Philippines')
            x_min = 110;
            x_max = 180;
            y_min =  -20.0;
            y_max = 35.0;
            
        elseif (strcmp(reg_lab,'nwpac'))
            disp('northwest pacific near taiwan')
            y_min = 14;
            y_max = 26;
            x_min = 117;
            x_max = 140;
            
        elseif (strcmp(reg_lab,'ccs'))
            disp('california current region')
            y_min = 25;
            y_max = 50;
            x_min = 210;
            x_max = 255;
            
        elseif (strcmp(reg_lab,'bob'))
            disp('bay of bengal')
            y_min = -5;
            y_max = 24;
            x_min = 77;
            x_max = 101;
            
        elseif (strcmp(reg_lab,'gom'))
            disp('gulf of mexico')
            y_min = 10;
            y_max = 30;
            x_min = 260;
            x_max = 310;
            
        elseif(strcmp(reg_lab,'red'))
            disp('red sea')
            y_min = 7;
            y_max = 32;
            x_min = 31;
            x_max = 78;
            
        elseif(strcmp(reg_lab,'ioc'))
            disp('ioc')
            y_min = -28;
            y_max = 32;
            x_min = 28;
            x_max = 122;
            
        elseif(strcmp(reg_lab,'hawaii'))
            disp('hawaii')
            y_min = 20;
            y_max = 40;
            x_min = 180;
            x_max = 240;
            
        elseif(strcmp(reg_lab,'SO_31S'))
            disp('SO_31S')
            y_min = -31;
            y_max = -31;
            x_min = 0;
            x_max = 360;
            
        elseif (strcmp(reg_lab,'bigeqpac'))
            disp('whole tropical pacific')
            x_min = 116;
            x_max = 290;
            y_min = -13.0;
            y_max =  13.0;
            
        elseif (strcmp(reg_lab,'ceckr'))
            disp('ceckr')
            x_min = 125;
            x_max = 135;
            y_min = 30.0;
            y_max =  40.0;
            
        end
        
        disp(['select in lat region ' num2str(y_min) ...
            ' to ' num2str(y_max)])
        disp(['select in lon region ' num2str(x_min) ...
            ' to ' num2str(x_max)])
        
        %% tolerance for region
        %% makes sure discretization doesn't mess us up
        %% this should be about a grid point in size
        %% use two, just in case...
        %% also makes it easier to compare with old files
        x_eps = 0.50;
        y_eps = 0.50;
        
        %% see if there is subsampling
        if (max(i_sub_xy)>1)
            disp(['subsample factor in x,y: ' num2str(i_sub_xy)])
            disp(['caution!  This may cut off top and left edge of region'])
            disp('increase tolerance by subsample factor')
            if (i_sub_xy(1) > 1)
                x_eps = x_eps*i_sub_xy(1);
            end
            if (i_sub_xy(2) > 1)
                y_eps = y_eps*i_sub_xy(2);
            end
        end
        
        %% looping over all year from 1993 - 2016
        years = [1993:2020];lyr = length(years);
        mon = [1:12]; lmn = length(mon);
        
        %% increment
        cc = 0;
        
        for yy = 1:lyr
            y1 = years(yy);
            dpath1 = [dpath num2str(y1) '/'];
            
            if y1 == 2020
                mon = [1:5]; lmn = length(mon);
            end
            
            for mm = 1:lmn
                m1 = mon(mm);
                if m1 < 10
                    mlab = ['0' num2str(m1)];
                else
                    mlab  = num2str(m1);
                end
                dpath2 = [dpath1 mlab '/'];
                %% choose base filename for aviso
                %% depends on grid type
                %% choose whether to use quarter degree or not
                disp('read 1/4 degree lat/long grid files')
                namfile = 'dt_global_allsat_phy_l4_';
                
                %% number of characters before date strings
                nchar = 32;
                
                %% read in a list of file names (on his disk)
                %% becomes a cell array
                files = importdata([dpath2 'dataset-duacs-rep-global-merged-allsat-phy-l4-v3_' num2str(y1) '_' mlab '.txt']);
                
                %% convert the input file list
                %% cells to characters
                filenames = char(files);
                %% Remove whitespaces
                filenames = strtrim(filenames);
                
                %% --> Set model/data timing.
                %% could just get all files that are there, but the following code
                %% lets you choose the first and last day to use.
                
                %% choose first date to get
                
                %% delayed time data
                date1 = datenum(y1,m1,01,0,0,0);
                date2 = datenum(y1,m1,eomday(y1,m1),0,0,0);
                
                %% CAUTION:: 2019 DT is only upto 20191015
                %% if y1 == years(end) & m1 == 10
                %%     date2 = datenum(y1,10,15,0,0,0);
                %% end
                
                %% have data at 1 day intervals
                deltat = 1;
                
                %% now date1 is the first day to get, date2 is the last
                %% make Data files dates extension.
                %% matlab datenum (used for output of .mat file)
                tsave1 = [date1:deltat:date2];
                
                %% date string for filename
                tabdat = datestr(tsave1',30);
                
                %% cut off the time part, so only 8 characters
                tabdat = tabdat(:,1:8);
                
                %% Compute number of days.
                %% should be first dimension of tabdat
                ndys = size(tabdat,1);
                disp([num2str(ndys) ' days']);
                
                %% choose whether or not to plot daily data
                %% iplot_day = 1;
                iplot_day = 0;
                
                %% --> Read data.
                for nt = 1:ndys
                    %% nt = 1;
                    %% Read aviso data field.
                    %% look for file in list; compare first nchar characters
                    %% (so don't have to know the date when updated)
                    %% have to do it on cell array!
                    
                    nfile = find(strncmp(files,[namfile tabdat(nt,:)],nchar)==1);
                    
                    if (length(nfile) < 1)
                        disp(['nfile < 1 in aviso_glb_sla.m (file missing??)'])
                        keyboard
                    elseif (length(nfile) > 1)
                        disp(['nfile > 1 in aviso_glb_sla.m (file duplicated??)'])
                        disp(filenames(nfile,:))
                        disp('take the later one and continue')
                        nfile = nfile(end);
                    end
                    
                    disp([' reading ' filenames(nfile,:)]);
                    
                    %% using ncload
                    %% ncload([dpath filenames(nfile,:)]);
                    %% dont need to use old ncload.m
                    
                    %% cpapadop
                    %% using new matlab netcdf internals
                    ncid = netcdf.open([dpath2 filenames(nfile,:)],'NC_NOWRITE');
                    
                    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
                    for varid =1:numvars
                        varlist(varid) = {netcdf.inqVar(ncid,varid-1)};
                    end
                    
                    prec='single';
                    param = 'latitude';
                    ind=find(ismember(varlist,param));
                    yd = netcdf.getVar(ncid,ind-1,prec);
                    
                    param = 'longitude';
                    ind=find(ismember(varlist,param));
                    xd = netcdf.getVar(ncid,ind-1,prec);
                    
                    %% first dimension is latitude
                    %% sea surface height in
                    scale_factor = 10^-4;
                    if madt
                        param = 'adt';
                        ind=find(ismember(varlist,param));
                        temp_ssh = netcdf.getVar(ncid,ind-1,prec);
                        
                        fillval = min(temp_ssh(:));
                        temp_ssh(temp_ssh == fillval) = nan;
                        
                        ssh = (temp_ssh.*scale_factor)';
                    elseif msla
                        param = 'sla';
                        ind=find(ismember(varlist,param));
                        temp_ssh = netcdf.getVar(ncid,ind-1,prec);
                        
                        fillval = min(temp_ssh(:));
                        temp_ssh(temp_ssh == fillval) = nan;
                        
                        ssh = (temp_ssh.*scale_factor)';
                    end
                    netcdf.close(ncid);
                    
                    if ( cc == 0 && nt == 1)
                        disp('first data file')
                        %% Select data in the model domain
                        %% select points in desired region (range)
                        ix = find(xd >= (x_min-x_eps) & xd <= (x_max+x_eps));
                        iy = find(yd >= (y_min-y_eps) & yd <= (y_max+y_eps));
                        
                        %% for 31S for Lynn, we need no tolerance
                        if(strcmp(reg_lab,'SO_31S'))
                            iy = min(find(yd >= y_min));
                        end
                        
                        %% see if there is subsampling
                        if (max(i_sub_xy)>1)
                            if (i_sub_xy(1) > 1)
                                ix = ix(1:i_sub_xy(1):end);
                            end
                            if (i_sub_xy(2) > 1)
                                iy = iy(1:i_sub_xy(2):end);
                            end
                        end
                        
                        figure
                        imagesc(xd(ix),yd(iy),ssh(iy,ix))
                        axis xy
                        colorbar
                        title('first ssh field')
                        
                        %% create arrays
                        nx = length(ix);
                        ny = length(iy);
                        %% file to save data (dimensioned: time, y, x)
                        yr_ln = lyr.*365;
                        dsave = nan(yr_ln,ny,nx);
                        tsave = nan(yr_ln,1);
                        %% reducd axes
                        xsave = xd(ix);
                        ysave = yd(iy);
                    end
                    if (iplot_day == 1)
                        figure(2)
                        imagesc(xsave,ysave,ssh(iy,ix))
                        axis xy
                        colorbar
                        title([num2str(nt) '-th day ssh field'])
                    end
                    %% put ssh in save array (note: y is first)
                    cc = cc + 1;
                    dsave(cc,:,:) = ssh(iy,ix);
                    tsave(cc,1) = tsave1(1,nt);
                end
            end
        end
        clear clear lat lon adt
        
        %% cleaning dsave
        in = find(~isnan(tsave));
        dsave = squeeze(dsave(in,:,:));
        tsave = squeeze(tsave(in,1));
        
        %% choose whether or not to save to a mat file
        isave = 1;
        
        if (isave == 1 && madt == 1)
            %% save aviso dsave xsave ysave tsave
            matfname = [opath 'aviso_tot_MADT_' reg_lab];
            
        elseif (isave == 1 && msla == 1)
            matfname = [opath 'aviso_tot_MSLA_' reg_lab];
            
        end
        
        disp(['save to mat file ' matfname])
        %% save([matfname],'dsave','xsave','ysave','tsave')
        
        %% for big files (foe ex: trop pac)
        save([matfname],'-v7.3','dsave','xsave','ysave','tsave')
    end
end
diary off;
return

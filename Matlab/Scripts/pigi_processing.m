function [out_data,comnt,comnt_t] = pigi_processing(data_dir,dep);

%--------------------------------------------------------------------------
% Function to read in optode/GTD data PIGI output LabVIEW files. Compiles
% all data from a single acquisition.
% NOTE: DOES NOT perform calculations, or post-processing on data 
% 
% INPUT:
% data_dir = directory in which data is saved
% dep = structure containing deployment info (e.g. dates, location, etc.)
%     dep..name = deployment name
%     dep.year = year of deployment
%     dep.start = start of deployment (datestring format)
%     depo.finish = end of deployment (datestring format)
%
% OUTPUT:
% out_data = data structure
%     out_data.time = Matlab time format
%     out_data.gps_time = Matlab time format
%     out_data.gps_lat = degrees N (if GPS connected to PIGI)
%     out_data.gps_long = degrees E (if GPS connected to PIGI)
%     %Raw Optode / GTD data
%         out_data.raw.o2uM = O2 concentration (umol/L)
%         out_data.raw.o2sat = O2 saturation (%)
%         out_data.raw.tp = GTD total dissolve gas presure (mbar)
%         %Raw instrument signals - optode (refer to manual for details):
%             out_data.raw.opt_calPhase = cal_phase
%             out_data.raw.opt_tcPhase = TC_phase
%             out_data.raw.opt_c1rph = c1rph
%             out_data.raw.opt_c2rph = c2rph
%             out_data.raw.opt_c1amp = c1_amp
%             out_data.raw.opt_c2amp = c2_amp
%             out_data.raw.opt_rawT_V = opt_raw_T
%         %Raw instrument signals - GTD (refer to manual for details):
%             out_data.raw.gtd_supplyV = supply_V
%             out_data.raw.gtd_ana1 = ana1
%             out_data.raw.gtd_ana2 = ana2
%             out_data.raw.gtd_dig1 = dig1
%             out_data.raw.gtd_dig2 = dig2        
%     out_data.opt_T = Optode temperature (C)
%     out_data.gtd_T = GTD temperature (C)
%     out_data.inst_flow = Instruments flow rate (L/min)
%     out_data.prim_flow = Primary chamber flow rate (L/min)
%     out_data.flow_status = Pump status (1 = on)%     
%     out_data.units = units for all variables
% comnt = string array of underway comments recorded during the acquisition
% comnt_t = array of time stamps corresponding with each underway comment
% 
% R. Izett
% rizett@eoas.ubc.ca
% UBC Oceanography
% Last Updated: June 2020
%--------------------------------------------------------------------------

cd(data_dir);

display('Extracting data. Please wait...')

%--- Specify number of lines of each new data file to remove at the start 
    %(e.g. remove bad data at the beginning of an acquisition / deployment)
    remove = 10;
    
%--- GET FILENAMES
    % create list of filenames
        filz = dir('*.lvm');
        for ff = 1:length(filz); fnames{ff,:} = filz(ff).name; end; clear ff filz
    
    %number of files
        [r,c]=size(fnames);
        
%--- Create list of headers and initialize data variables
    %Potential headers:
        pot_headerz = {'comp_jd' 'o2sat' 'o2' 'o2_T' 'pres' 'pres_T' 'flow_instruments' 'gps_jd' 'gps_lat' 'gps_long' ...
            'cal_phase' 'TC_phase' 'c1rph' 'c2rph' 'c1_amp' 'c2_amp' 'opt_raw_T' 'supply_V' 'ana1' 'ana2' 'dig1' 'dig2' ...
            'pump_status' 'flo2' 'pres2' 'pres_T2' 'prime_status' 'prime_status'};
        
    % Load headers names from file
        fid = fopen(['PIGI_headers.txt'],'rt'); %open .dat file
        headerz = textscan(fid,repmat('%s',1,1),100,'headerlines',0,'collectoutput',1); headerz = headerz{1};
        fclose(fid);
        
    %See if potential headers exist
        for kk = 1:numel(pot_headerz)
            is_here = strcmp(pot_headerz{kk}, headerz);
            
            if all(is_here == 0); %if the header doesn't exist, make the data nan
                s = sprintf('%s = ''variable not present'';',pot_headerz{kk}); eval(s); clear s
            end
            clear is_here 
        end     
        
    %Initialize variables to hold data and underway comments
        colm = length(headerz);
        dat = nan(1,colm);
        comnt = []; %comments 
        comnt_t = []; %comments time
 
%--- READ DATA AND CONCATENATE   
    for ii=1:r; %go through each file
        
        fi = dir([fnames{ii,:}]);
        if fi.bytes < 4500; continue; end
        fid = fopen(char(fnames(ii,:)),'rt'); %open .dat file
        readdat = textscan(fid,[repmat('%f',1,colm),repmat('%s',1,100)],'headerlines',23,'collectoutput',1); catdat = readdat{1}; stringz = readdat{2}; clear readdat
        fclose(fid);

        %get comments from file
        for kk = 1:length(stringz(:,1))
            if ~isempty(stringz{kk,1})
                temp = '';
                    for jj = 1:100
                        temp = strcat(temp,stringz(kk,jj), string(' '));
                    end
                comnt = [comnt; temp];
                comnt_t = [comnt_t,catdat(kk,1)];
                clear temp
            end
        end
        clear stringz kk        
        
        %remove first X rows from data
            if remove > length(catdat(:,1))
                clear fid catdat
                continue;
            else
                catdat(1:remove,:)=[];
            end
        
        %concatonate data
        if ~isempty(catdat)
            dat=cat(1,dat,catdat); 
            dat=[dat;nan(1,length(catdat(1,:)))];
            
        end
        clear fid catdat
    end
   
    dat(1,:)=[];  %remove first row of NaNs used to initialize the matrix
        
%--- SORT DATA BY TIME
    [T,I] = sort(dat(:,1)); clear T
    dat = dat(I,:);
    
%--- CREATE NAMED VARIABLES
    for hh = 1:length(headerz)
        s = sprintf('%s = dat(:,hh);',headerz{hh,:}); eval(s);
    end; clear hh 

%--- DO A QUICK CLEANUP (DE-SPIKE) ON DATA
    o2sat(o2sat > 200) = nan;
    o2sat(o2sat < 20) = nan;
    o2_T(o2_T > 40) = nan;
    o2(o2 < 50) = nan;
    o2(o2 > 700) = nan;
    pres(pres < 500) = nan;
    pres(pres > 1600) = nan;
  
%--- CONVERT JULIAN DAYS TO MATLAB TIME
    mdate = datenum(dep.year, 1, comp_jd); 
    gps_mdate = datenum(dep.year, 1, gps_jd); 
    comnt_t = datenum(dep.year, 1, comnt_t);
    
%--- Replace missing data w/ nans
     for kk = 1:numel(pot_headerz)
        is_here = strcmp(pot_headerz{kk}, headerz);
            
        if all(is_here == 0); %if the header doesn't exist, make the data nan
            s = sprintf('%s = nan(size(mdate));',pot_headerz{kk}); eval(s); clear s
        end
        clear is_here 
    end     
    
%--- make data structure
    out_data.time = mdate;
    out_data.gps_time = gps_mdate;
    out_data.gps_lat = gps_lat;
    out_data.gps_long = gps_long;
    
    %Raw Optode / GTD data
        out_data.raw.o2uM = o2; %concentration, uM
        out_data.raw.o2sat = o2sat; %percent-saturation, %/100
        out_data.raw.tp = pres; %dissolved gas pressure, mbar
        %Raw instrument signals - optode (refer to manual for details)
        out_data.raw.opt_calPhase = cal_phase;
        out_data.raw.opt_tcPhase = TC_phase;
        out_data.raw.opt_c1rph = c1rph;
        out_data.raw.opt_c2rph = c2rph;
        out_data.raw.opt_c1amp = c1_amp;
        out_data.raw.opt_c2amp = c2_amp;
        out_data.raw.opt_rawT_V = opt_raw_T;
        %Raw instrument signals - GTD (refer to manual for details)
        out_data.raw.gtd_supplyV = supply_V;
        out_data.raw.gtd_ana1 = ana1;
        out_data.raw.gtd_ana2 = ana2;
        out_data.raw.gtd_dig1 = dig1;
        out_data.raw.gtd_dig2 = dig2;
        
    out_data.opt_T = o2_T;
    out_data.gtd_T = pres_T;
    out_data.inst_flow = flow_instruments; %flow over instruments
    out_data.prim_flow = flo2; %flow through primary loop
    out_data.flow_status = pump_status; %pump on or off
    
    out_data.units.time = 'Matlab format';
    out_data.units.o2uM = 'umol/L';
    out_data.units.o2sat = '%/100';
    out_data.units.tp = 'mbar';
    out_data.units.T = 'deg-C';
    out_data.units.flow = 'L/min';
   
    display(' ')
    display('Data extracted');
    display(' ')
    
end


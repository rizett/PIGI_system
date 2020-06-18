function [sname] = save_version(sdir,sname,data,varname)

%--------------------------------------------------------------------------
% Function to save data to a specified directory (dir) with a version / 
% date appended (i.e. name_vYYYYMMDD.mat)
%
% USAGE:
%   save_version(sdir,sname,data,varname)
% 
% INPUTS:
%   sdir = directory where data will be saved
%   sname = file name (without extension)
%   data = the variable to be saved
%   varname = variable name of data to be saved 
% 
% OUTPUTS:
%   sname = directory location of saved data
% 
% Last updated: June 2020
% R. Izett, rizett@eoas.ubc.ca
% UBC Oceanography
%--------------------------------------------------------------------------

%--- create a new variable name for saving
    save_name = [sname,'_v',datestr(datenum(date),'yyyymmdd'),'.mat'];

%--- create save directory / name    
    if ispc
        sname = [sdir,'\',save_name]; %PC
    else
        sname = [sdir,'/',save_name]; %Mac
    end

%--- save data
    varname_string = varname;
    eval(sprintf('%s = data;',varname));

    save(sname,varname_string);

return
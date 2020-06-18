function [c_o2,s_o2,p_o2] = optode_recalc(calphase,new_T,S,slp,fname,func)

%-------------------------------------------------------------------------
% Recalculate optode cO2, sO2, pO2 using specified function (svu =
% Stern-Volmer Uchida; ssv = Simplified Stern-Volmer; poly = emperical
% polynomial). Salinity compensation is also performed.
% 
% INPUT:
%   calphase = raw Optode phase shift (calPhase - see Optode manual)
%   new_T = Sea surface temperature (deg-C)  used to re-calculate O2 (best 
%       to use TSG or CTD-corrected underway SST) 
%   S = Salinity (PSU); for S-compensation and unit conversions 
%   slp = sea level pressure (mbar)
%   fname = name of file containing optode saved settings 
%       (.txt file containing output from running "get all\r\n" on Optode)
%       use full file directory (e.g. /file/location/optode_settings.txt)
%   func = desired function for calculating O2 (default = svu)
%       'svu' = Stern-Volmer Uchida (Uchida et al., 2008)
%           (*NOTE: requires that Optode foil coefficients have been 
%           obtained using a multi-point calibration)
%       'ssv' = Simplified Stern-Volmer (GEOMAR/Bittig et al.,2018); 
%           (*NOTE: requires that Optode fiol coefficients have been derivd 
%           using this formula during mulit-point  calibration)
%       'poly' = emperical polynomial 
%           (Use this function if no foil coefficients are provided / if
%           the sensor has not been multi-point calibrated)
% 
% OUTPUT:
%   c_o2, s_o2, p_o2 = re-calculated O2 concentration [uM], O2 saturation
%   [%/100], and partial pressure [mbar]
% 
% Last updated: June 2020
% R. Izett, rizett@eoas.ubc.ca
% UBC Oceanography
%--------------------------------------------------------------------------

%--- Set default function
    if nargin < 6
        warning('Function not specified. ''svu'' set by default')
        func = 'svu';
    end
    
%--- Set default SLP
    if ~exist('slp','var')
        slp = 1013.25;
    elseif isempty(slp);
        slp = 1013.25;
    end

%--- Scan Optode settings from file
    fid = fopen(fname,'rt');
    settings = textscan(fid,repmat('%s',1,1),200,'headerlines',0,'collectoutput',1); 
    settings = settings{1};
    fclose(fid);

%---Coefficients for salinity compensation 
    %See Optode manual
    B0 = -7.01577e-3;
    B1 = -7.70028e-3;
    B2 = -1.13864e-2;
    B3 = -9.51519e-3;
    C0 = -2.75915e-7;
    
    ts = log((298.15-new_T)./(273.15+new_T));

%--- Re-calculate O2:
    %SVU equation
    %Uchida et al., 2008. See Optode manual.
        if strcmp(func,'svu')
            %Get coefficients from file
                sv = settings([153:159]);
                ksv = str2num(sv{1}) + str2num(sv{2}) .* new_T + str2num(sv{3}) .* new_T.^2;
                po = str2num(sv{4}) + str2num(sv{5}) .* new_T;
                pc = str2num(sv{6}) + str2num(sv{7}) .* calphase;

            %Calculate O2 conc.
                c_o2 = (po./pc - 1)./ksv; %uM

        end

    %Simplified SV
        if strcmp(func,'ssv')
            %Get coefficients from file
                sv = settings([153:159]);
                ksv = str2num(sv{1}) + str2num(sv{2}) .* new_T + str2num(sv{3}) .* new_T.^2;
                po = 1 + str2num(sv{4}) .* new_T;
                pc = str2num(sv{5}) + str2num(sv{6}) .* calphase + str2num(sv{7}) .* calphase.^2;

            %Calculate pO2 and O2 conc.
                p_o2 = (po./pc - 1)./ksv; %mbar
                c_o2 = O2ptoO2c(p_o2,new_T,S,0);           
           
        end

    %Empirical / 27th order polynomial
        if strcmp(func,'poly')
            %Get coefficients from file
                c = settings([57:70,74:87]); %C0-C13 = FoilCoefA; C14-C27 = FoilCoefB
                m = settings(91:118); %m0-m27 = FoilPolyDegT
                n = settings(122:149); %n0-n27 = FoilPolyDeg0
            
            %Calculate pO2 and O2 conc.
                p_o2 = 0;
                for kk = 1:numel(c)
                    p_o2 = p_o2 + (str2num(c{kk})) .* (new_T .^ str2num(m{kk})) .* (calphase .^ str2num(n{kk}));
                end

                c_o2 = O2ptoO2c(p_o2,new_T,S,0);

        end
        
%--- Salinity compensation
    %salinity component of O2 solubilty
        scorr = exp(S.*(B0 + ts .* B1 + B2 .* ts.^2 + B3 .* ts.^3) + C0.*S.^2); 
    % saturation water vapor in mbar, 0 PSU
        ph2osat_0 = 1013.25 .* (exp(24.4543 - (67.4509 .* (100 ./ (new_T+273.15))) - (4.8489 .* log(((273.15 + new_T) ./ 100))) -0.000544 .* repmat(0,size(new_T)))); 
    % saturation water vapor in mbar
        ph2osat = 1013.25 .* (exp(24.4543 - (67.4509 .* (100 ./ (new_T+273.15))) - (4.8489 .* log(((273.15 + new_T) ./ 100))) -0.000544 .* S)); 
    %salinity dependence of vapour pressure
        vcorr = (1013.25 - ph2osat_0) ./ (1013.25 - ph2osat); 
    %full Salinity compensation on O2 (as in Bittig et al., 2018)
        c_o2 = c_o2 .* scorr .* vcorr; 

    %Unit conversions
        s_o2 = O2ctoO2s(c_o2,new_T,S,0,slp); %saturation (%/100)
        p_o2 = O2ctoO2p(c_o2,new_T,S,0); %partial pressure (mbar)   
        
return

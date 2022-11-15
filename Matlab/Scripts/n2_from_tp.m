function [dat,p,air] = n2_from_tp(tp,o2,ar,pco2,T,S,air)

%--------------------------------------------------------------------------
%Calculate N2 from TP and other available gas data.
%
% INPUT:
%     tp = total dissolved gas pressure (temperature and bias corrected); mbar
%     o2 = calibrate O2 saturation; %/100
%     ar = Argon saturation; %/100 (leave empty, [], if not available)
%     pco2 = Partial pressure of CO2; ppm (leave empty, [], if not available)
%     T = SST; C
%     S = salinity; PSU
%     air = data structure containing:
%        air.slp = sea level pressure; mbar
%
% OUTPUT:
%     dat = data structure containing:
%         dat.n2sat = Nitrogen saturation; %/100
%         dat.n2_molkg = Nitrogen concentration; mol/kg
%         dat.do2n2 = delta-O2/N2; %
%         dat.units = units for all variables
%     p = data structure containing (partial) pressure data (all mbar):
%         p.o2 = Oxygen partial pressure
%         p.ar = Argon partial pressure
%         p.co2 = CO2 partial pressure
%         p.h2o = saturation vapour pressure
%     air = data structure containing:
%        air.slp = sea level pressure; mbar
%        air.ph2o = vapour pressure; mbar
%        air.x_o2 = atmospheric O2 mixing ratio
%        air.x_ar = atmospheric Ar mixing ratio
%        air.x_n2 = atmospheric N2 mixing ratio
%
% Last updated: June 2020
% R. Izett, rizett@eoas.ubc.ca
% UBC Oceanography
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE N2 FROM O2 & TP DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--- Set atmospheric mixing ratios
        air.x_o2 = 0.20946;
        air.x_n2 = 0.78084;
        air.x_ar = 0.0093;
        
    %--- Atmospheric pH2O (mbar)
        air.ph2o = svp(T); %Assume 100 % humidity at air-water interface;
        
    %--- Parial pressure calculations (all mbar)
        p.o2 = O2stoO2p(o2.*100,T,S,0,air.slp);

        %Argon
        if ~isempty(ar)
            p.ar = ar .* air.x_ar .* (air.slp - air.ph2o); 
        end
        
        %CO2, ppm
        if isempty(pco2)
            p.co2 = repmat(400,size(p.o2))*1e-6;
        else
            p.co2 = pco2 .* 1e-6 .* (air.slp-air.ph2o); %convert uatm to mbar
        end
        
        %Vapour pressure in seawater
        p.h2o = vpress(S,T).*1013.25; %Weiss & Price 1980 saturation vapour pressure         

    %--- Calculate pN2
        if isempty(ar)
            p.n2 = (tp - p.o2 - p.h2o - p.co2) / (1 + ((1 - air.x_n2 - air.x_o2)/air.x_n2));
        else
            p.n2 = (tp - p.o2 - p.h2o - p.ar - p.co2) / (1 + ((1 - air.x_n2 - air.x_o2 - air.x_ar)/air.x_n2));
        end
                
    %--- N2 saturation level [%]
        dat.n2sat = p.n2 ./ (air.x_n2 .* (air.slp - air.ph2o)); 
        
    %--- N2 concentration [mol/kg]
        dat.n2_molkg = dat.n2sat .* N2sol(S,T) ./1e6 .* (air.slp - air.ph2o) ./ (1013.25 - air.ph2o);
       
    %--- delO2/N2
        dat.do2n2 = (o2 ./ dat.n2sat - 1)*100;
        
    dat.units.n2sat = '%/100';
    dat.units.n2_molkg = 'mol/kg';
    dat.units.do2n2 = '%';
    
    p.units = 'mbar';
    
return        
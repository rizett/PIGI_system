function SVP = svp(T)

%--------------------------------------------------------------------------
% Calculate saturation vapour pressure from atmospheric temperature.

% T = atmospheric temperature [C]
% SVP = saturation vapor pressure according to Goff & Gratch Equation (1946)
% SVP in mbar
%--------------------------------------------------------------------------

T = T + 273.15;

Ts = 373.15;    % standard temperature at steam point on Kelvin scale (Goff, 1965)
Ps = 101324.6;  % standard atmospheric pressure at steam point (Pascal)

log10SVP = -7.90298*(Ts./T - 1) + 5.02808*log10(Ts./T) -1.3816e-7*(10.^(11.344*(1 - T./Ts)) - 1) + 8.1328e-3*(10.^(3.49149*(1 - Ts./T)) - 1) + log10(Ps);

SVP = 10.^(log10SVP);

SVP = SVP * 0.01; %convert Pa --> mbar
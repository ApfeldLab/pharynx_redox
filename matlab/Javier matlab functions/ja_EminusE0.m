function [EminusE0] = ja_E0(OxD,E,z)

if nargin <3
    z=2;
end
if nargin <2
    E= 0;
end

Temp = 22;
EminusE0= - (8314.462*(273.15+Temp)./(z.*96485.3415))*log((1-OxD)./OxD);
EminusE0 = real(EminusE0);
end




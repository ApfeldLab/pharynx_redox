function [E] = ja_E(OxD,E0,z,Temp)
% [E] = ja_E(OxD,E0,z,Temp)
% updated June 18 2013

if nargin <4
   Temp = 22;
end

if nargin <3 | z==[]
    z = 2;
end

if nargin <2 | E0==[]
    E0 = -265;
end

if sum([(max(size(OxD))) > 1 (max(size(E0))) > 1 (max(size(z))) > 1 (max(size(Temp))) > 1])> 1
    error('more than one multidimensional parameter')
else
      
    E = E0 - (8314.462*(273.15+Temp)./(z.*96485.3415))*log((1-OxD)./OxD);
    E = real(E);
end




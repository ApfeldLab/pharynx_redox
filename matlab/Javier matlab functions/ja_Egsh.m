function [Egsh] = ja_Egsh(OxD,GSHtotal,E0,z,Temp)
% [Egsh] = ja_Egsh(OxD,GSHtotal,E0,z,Temp)
% updated August 16 2013

if nargin <5
   Temp = 22;
end

if nargin <4 | z==[]
    z = 2;
end

if nargin <3 | E0==[]
    E0 = -240;
end

if nargin <2  | GSHtotal==[]
    GSHtotal = 0.01;
end

if sum([(max(size(OxD))) > 1 (max(size(GSHtotal))) > 1 (max(size(E0))) > 1 (max(size(z))) > 1 (max(size(Temp))) > 1])> 1
    error('more than one multidimensional parameter')
else
      
    Egsh = E0 - (8314.462*(273.15+Temp)./(z.*96485.3415))*log(2.*GSHtotal.*((1-OxD).^2)./OxD);
    Egsh = real(Egsh);
end




function [testoxd] = ja_E2oxd_gsh_withGSHtotal(E,GSHtotal,E0,z,Temp)
% [testoxd] = ja_E2oxd_gsh_withGSHtotal(E,GSHtotal,E0,z,Temp)
% This function calculates OxD_gsh as a function of E and E0
% ** Unlike ja_E2oxd_gsh, the effect of [GSH]total is NOT folded into E0 **
% Instead: E0 = -240 - (RT/zF) ln([GSH]total)
% Temp is in Celcius
% updated December 2 2013

if nargin <5
   Temp = 22;
end

if nargin <4 | isempty(z)
    z = 2;
end

if nargin <3 | isempty(E0)
    E0 = -240; %midpoint potential of GSSG:GSH couple
end

if nargin <2  | isempty(GSHtotal)
    GSHtotal = 0.01;
end

if sum([(max(size(E))) > 1 (max(size(GSHtotal))) > 1 (max(size(E0))) > 1 (max(size(z))) > 1 (max(size(Temp))) > 1])> 1
    error('more than one multidimensional parameter')
else
    
r = 8314.462;
t = 273.15+Temp;
f = 96485.3415;

a= r.*t./z./f;
E0 = E0 - a.* log(GSHtotal);

b = exp(-(E-E0)./a)/2;

testoxd = 0.5* ((2+b) - sqrt(b.*(b+4)));
testoxd = real(testoxd);

end

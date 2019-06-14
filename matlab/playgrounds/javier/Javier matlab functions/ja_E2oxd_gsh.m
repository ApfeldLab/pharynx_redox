function [testoxd] = ja_E2oxd_gsh(E,E0,z,Temp)
% [testoxd] = ja_E2oxd_gsh(E,E0,z,Temp)
% This function calculates OxD_gsh as a function of E and E0
% ** Note that the effect of [GSH]total is folded into E0 **
% E0 = E0' - (RT/zF) ln([GSH]total)
% Temp is in Celcius
% updated August 22 2013

if nargin <4
   Temp = 22;
end

if nargin <3 | z==[]
    z = 2;
end

if nargin <2 | E0==[]
    E0 = -240; %midpoint potential of GSSG:GSH couple
end

if sum([(max(size(E))) > 1 (max(size(E0))) > 1 (max(size(z))) > 1 (max(size(Temp))) > 1])> 1
    error('more than one multidimensional parameter')
else
r = 8314.462;
t = 273.15+Temp;
f = 96485.3415;


a= r.*t./z./f;

b = exp(-(E-E0)./a)/2;

testoxd = 0.5* ((2+b) - sqrt(b.*(b+4)));

end

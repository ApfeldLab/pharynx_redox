function [testoxd] = ja_E2oxd(E,E0,z,Temp)
% [testoxd] = ja_E2oxd(E,E0,z,Temp)
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

if sum([(max(size(E))) > 1 (max(size(E0))) > 1 (max(size(z))) > 1 (max(size(Temp))) > 1])> 1
    error('more than one multidimensional parameter')
else
r = 8314.462;
t = 273.15+Temp;
f = 96485.3415;


a= r.*t./z./f;

testoxd = (1./(1+exp(-(E-E0)./a)));

end

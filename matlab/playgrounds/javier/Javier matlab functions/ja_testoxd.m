function [testoxd] = ja_testoxd(E,E0,z)

if nargin <3
    z=2;
end
if nargin <2
    E0= -265;
end

Temp = 22;
r = 8314.462;
t = 273.15+Temp;
f =96485.3415;


a= r*t/z/f;

testoxd = (1./(1+exp(-(E-E0)./a)));

end

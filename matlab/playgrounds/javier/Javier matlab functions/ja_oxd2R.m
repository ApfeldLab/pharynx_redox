function [testR] = ja_oxd2R(oxd,Rmin,Rmax,Instrument_Factor)
% [testR] = ja_oxd2R(oxd,Rmin,Rmax,Instrument_Factor)
% updated June 18 2013

if nargin <4
    Instrument_Factor = 0.171;
end

if nargin <3 || Rmax==[]
   Rmax = 5.207;
end

if nargin <2 || Rmin==[]
    Rmin = 0.667;
end


if sum([(max(size(oxd))) > 1 (max(size(Rmin))) > 1 (max(size(Rmax))) > 1 (max(size(Instrument_Factor))) > 1])> 1
    error('more than one multidimensional parameter')
else
    
testR = (Rmin + oxd .* (Instrument_Factor .* Rmax - Rmin)) ./ (1 - oxd .* (1- Instrument_Factor));

end

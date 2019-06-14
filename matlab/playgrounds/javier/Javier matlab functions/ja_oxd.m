function [oxd] = ja_oxd(R,Rmin,Rmax,Instrument_Factor)
% [oxd] = ja_oxd(R,Rmin,Rmax,Instrument_Factor)
% updated June 18 2013

if nargin <4
    Instrument_Factor = 0.171;
end

if nargin <3 | Rmax==[]
   %Rmax = 5.207
   Rmax = 6.65;
end

if nargin <2 | Rmin==[]
    %Rmin = 0.667;
%      Rmin = 0.667*Rmax/5.207;
    Rmin = .852;
end

if sum([(max(size(R))) > 1 (max(size(Rmin))) > 1 (max(size(Rmax))) > 1 (max(size(Instrument_Factor))) > 1])> 1
    error('more than one multidimensional parameter')
else
    
oxd = (R-Rmin) ./((R-Rmin)+Instrument_Factor *(Rmax-R));

end




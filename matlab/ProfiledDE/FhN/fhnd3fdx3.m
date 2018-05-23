function r = fhnd2fdx3(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnd2fdy3
%
% The third derivative of the FitzHugh Nagumo equations with respect 
% to the input x. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell([2,2,2,3]);
r(1:2,1:2,1:2,1:2) = {0};
r{1,1,1,1} = -2*p(3);
end
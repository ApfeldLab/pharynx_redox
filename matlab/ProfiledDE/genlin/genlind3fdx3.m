function r = genlind3fdx3(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlind3fdx3
%
% The third derivative of linear differential equations with respect 
% to the input x. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(fd_cell,2);
r = cell(n,n,n,n);
r(:) = {0};

end
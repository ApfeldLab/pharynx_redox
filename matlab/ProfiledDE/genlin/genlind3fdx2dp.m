function r = genlind3fdx2dp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlind3fdx2dp
%
% The third derivative of linear differential equations with respect 
% to the input x (twice) and p(once) at time t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(fd_cell,2);
r = cell(n,n,n,length(p));
r(:) = {0};

end
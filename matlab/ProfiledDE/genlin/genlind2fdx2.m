function r = genlind2fdx2(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlind2fdx2
%
% The second derivative of linear differential equations with respect 
% to the input x at time t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = zeros(length(t),2,2,2);

n = size(fd_cell,2);
r = cell(n,n,n);
r(:) = {0};

end
function r = genlind3fdxdp2(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlind3fdxdp2
%
% The third derivative of linear differential equations with respect to 
% x (once) and p (twice) at times t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(fd_cell,2);
m = length(p);
r = cell(n,n,m,m);
r(:) = {0};

end
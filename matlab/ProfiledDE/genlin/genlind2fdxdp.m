function r = genlind2fdxdp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlind2fdxdp
%
% The second derivative of linear differential equations with respect to 
% x and p at times t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<4)
    more = [];
end

n = size(fd_cell,2);

more = checkmore(more,n);

r = cell(n,n,length(p));
r(:) = {0};

for i = 1:size(more.sub,1) 
    r{more.sub(i,1),more.sub(i,2),i} = 1;    
end

end
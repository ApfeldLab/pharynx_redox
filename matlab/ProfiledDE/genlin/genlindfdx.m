function r = genlindfdx(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlindfdx
%
% The derivative of linear differential equations with respect to x
% at times t with parameters p. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    more = [];
end

n = size(fd_cell,2);

more = checkmore(more,n);
pmat = more.mat;
pmat(sub2ind([n n],more.sub(:,1),more.sub(:,2))) = p(1:size(more.sub,1));

r = cell(n,n);

for i = 1:n
    for j = 1:n
        r{i,j} = pmat(i,j);
    end
end

end
function r = polyfun(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyfun
%
% One dimensional polynomial function with forcing -- vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);

r = cell(1);
r(:) = {0};

for i = 1:(length(p)-1) 
    r{1} = r{1} + p(i)*y{1}.^i;
end 

if nargin>=4 
    r{1} = r{1} + p(length(p))*more.forcing(t,more.fs);
else
    r{1} = r{1} + p(length(p));
end

end
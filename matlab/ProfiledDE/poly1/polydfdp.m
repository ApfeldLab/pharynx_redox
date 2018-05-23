function r = polydfdp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polydfdp
%
% One dimensional polynomial function with forcing
%
% Derivative with respect to parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);

r = cell(1,length(p));
r(:) = {0};

for i = 1:(length(p)-1)
    r{i} = y{1}.^i;
end

if nargin>=4 
    r{length(p)} = more.forcing(t,more.fs);
else
    r{length(p)} = 1;
end

end
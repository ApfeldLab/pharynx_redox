function r = polydfdx(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polydfdx
%
% One dimensional polynomial function with forcing
%
% Derivative with respect to component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);
r = cell(1);
r(:) = {0};

for i = 1:(length(p)-1)
    r{1} = r{1} + i*p(i)*y{1}.^(i-1);
end

end
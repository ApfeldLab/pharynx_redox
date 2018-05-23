function r = polyd3fdx3(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyd3fdx3
%
% One dimensional polynomial function with forcing
%
% Third derivative with respect to component. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(1,1,1,1);
r(:) = {0};

for i = 3:(length(p)-1) 
    r{1} = r{1} + i*(i-1)*(i-2)*p(i)*y{1}^(i-3);
end

end
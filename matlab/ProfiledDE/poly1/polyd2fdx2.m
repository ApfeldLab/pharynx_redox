function r = polyd2fdx2(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyd2fdx2
%
% One dimensional polynomial function with forcing
%
% Second derivative with respect to component. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);

r = y;

for i = 2:(length(p)-1) 
   r{1} = p(i)*i*(i-1)*y{1}.^(i-2); 
end

end
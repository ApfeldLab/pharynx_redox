function r = polyd3fdx2dp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyd3fdx2dp
%
% One dimensional polynomial function with forcing
%
% Third derivative with respect to input (twice) and parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y = eval_fdcell(t,fd_cell,0);
r = cell(1,1,1,length(p));
r(:) = {0};

for i = 2:(length(p)-1) 
   r{1,1,1,i} = i*(i-1)*y{1}.^(i-2);
end

end
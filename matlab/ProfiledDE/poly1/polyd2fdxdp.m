function r = polyd2fdxdp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyd2fdxdp
%
% One dimensional polynomial function with forcing
%
% Cross derivative with respect to input and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);

r = cell(1,1,length(p));
r(:) = {0};

for i = 1:(length(p)-1) 
    r{1,1,i} = i*y{1}.^(i-1);
end

end
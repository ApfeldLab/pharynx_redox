function r = forcingfun(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingfun
%
% forcing function estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


force_cell = Make_fdcell(p,more.basisp);

r = feval_check(more.fn,t,fd_cell,more.pars,more.extras);

fs = eval_fdcell(t,force_cell);

for i = 1:length(more.which)
    r{more.which(i)} = r{more.which(i)} + fs{i};
end



end
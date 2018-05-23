function r = forcingd3fdx3(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingdfdx
%
% forcing function estimation
%
% third derivative with respect to smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = feval_check(more.d3fdx3,t,fd_cell,more.pars,more.extras);


end
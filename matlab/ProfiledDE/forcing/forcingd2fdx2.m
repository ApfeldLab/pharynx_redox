function r = forcingd2fdx2(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingd2fdx2
%
% forcing function estimation
%
% second derivative with respect to smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = feval_check(more.d2fdx2,t,fd_cell,more.pars,more.extras);


end
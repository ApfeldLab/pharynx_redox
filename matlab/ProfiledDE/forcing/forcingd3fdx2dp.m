function r = forcingd3fdx2dp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingd3fdx2dp
%
% forcing function estimation
%
% third derivative with respect to smooth (twice) and forcing coefficents.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(size(fd_cell,2),size(fd_cell,2),size(fd_cell,2),length(p));
r(:) = {0};


end
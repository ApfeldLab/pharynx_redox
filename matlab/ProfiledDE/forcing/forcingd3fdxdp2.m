function r = forcingd3fdxdp2(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingd3fdxdp2
%
% forcing function estimation
%
% third derivative with respect to smooth and forcing coefficents (twice).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(size(fd_cell,2),size(fd_cell,2),length(p),length(p));
r(:) = {0};


end
function r = forcingd2fdxdp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingd2fdxdp
%
% forcing function estimation
%
% cross derivative with respect to smooth and forcing coefficents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(size(fd_cell,2),size(fd_cell,2),length(p));
r(:) = {0};


end
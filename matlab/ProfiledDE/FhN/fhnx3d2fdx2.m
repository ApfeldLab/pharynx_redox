function r = fhnx3d2fdx2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnx3d2fdx2
%
% Cubic approx to 1D FitzHugh Nagumo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(1,1,1);
r{1,1,1} = -2*p(2)*eval_fd(t,fd_cell{1});

end
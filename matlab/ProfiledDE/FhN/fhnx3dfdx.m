function r = fhnx3dfdx(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnx3dfdx
%
% Cubic approx to 1D FitzHugh Nagumo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = cell2mat(eval_fdcell(t,fd_cell,0)');
r = cell(1,1);

r{1,1} = p(1) - p(2)*eval_fd(t,fd_cell{1}).^2;

end
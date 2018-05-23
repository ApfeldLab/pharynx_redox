function r = lincubd2fdx2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnfun
%
% The FitzHugh-Nagumo equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(2,2,2);

y = eval_fdcell(t,fd_cell,0);
r(:) = {0};

r{1,1,1} = 6*p(7)*y{1};

end
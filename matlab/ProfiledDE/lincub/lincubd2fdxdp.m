function r = lincubd2fdxdp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnfun
%
% The FitzHugh-Nagumo equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(2,2,7);

y = eval_fdcell(t,fd_cell,0);
r(:) = {0};

r{1,1,1} = 1;
r{1,1,7} = 3*y{1}.^2;
r{1,2,2} = 1;
r{2,1,3} = 1;
r{2,2,4} = 1;

end
function r = lincubdfdp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnfun
%
% The FitzHugh-Nagumo equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);
r = cell(2,7);
r(:) = {0};

r{1,1} = y{1};
r{1,2} = y{2};
r{1,5} = 1;
r{1,7} = y{1}.^3;
r{2,3} = y{1};
r{2,4} = y{2};
r{2,6} = 1;

end
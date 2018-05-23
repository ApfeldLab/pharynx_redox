function r = lincubfdfdx(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnfun
%
% The FitzHugh-Nagumo equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell(2,2);

y = eval_fdcell(t,fd_cell,0);
r(:) = {0};

r{1,1} = p(1) + 3*p(7)*y{1}.^2;
r{1,2} = p(2);
r{2,1} = p(3);
r{2,2} = p(4);

end
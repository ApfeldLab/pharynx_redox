function r = lincubfun(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnfun
%
% The FitzHugh-Nagumo equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);
r = y;
r{1} = p(1)*y{1}+p(2)*y{2} + p(5) + p(7)*y{1}.^3;
r{2} = p(3)*y{1}+p(4)*y{2}+p(6);


end
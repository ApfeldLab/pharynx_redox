function r = heartfun(t,fd_cell,p)

y = eval_fdcell(t,fd_cell,0);
r = y;

r{1} = p(1) + p(2)*y{1} + p(3)*y{1}.^2 + p(4)*y{1}.^3 + p(5)*y{2} + p(6)*y{1}.*y{2};
r{2} = 0.06 + p(7)*y{1} + p(8)*y{2};

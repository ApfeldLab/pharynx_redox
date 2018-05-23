function r = heartdfdx(t,fd_cell,p)

y = eval_fdcell(t,fd_cell,0);
r = cell(2,2);

r{1,1} = p(2) + 2*p(3)*y{1} + 3*p(4)*y{1}.^2 + p(6)*y{2};
r{1,2} = p(5)+p(6).*y{1};

r{2,1} = p(7);
r{2,2} = p(8);


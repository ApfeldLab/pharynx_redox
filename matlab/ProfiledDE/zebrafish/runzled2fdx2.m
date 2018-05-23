function r = heartd2fdx2(t,fd_cell,p)

y = eval_fdcell(t,fd_cell,0);
r = cell(2,2,2);
r(:) = {0};

r{1,1,1} = 2*p(3) + 6*p(4)*y{1};
r{1,1,2} = p(6);
r{1,2,1} = p(6);


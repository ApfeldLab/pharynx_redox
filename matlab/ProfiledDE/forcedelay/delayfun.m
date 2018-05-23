function r = delayfun(t,fd_cell,p,force)


x = eval_fdcell(t,fd_cell,0);
f = force.beta(1) + force.beta(2)*monfn(t-p(3),force.force_fd);

r = x;

r{1} = p(1)*x{1} + p(2)*f;
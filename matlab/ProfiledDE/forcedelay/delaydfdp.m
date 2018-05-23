function r = delaydfdp(t,fd_cell,p,force)


x = eval_fdcell(t,fd_cell,0);
dforce = exp(eval_fd(t-p(3),force.force_fd));
force = force.beta(1) + force.beta(2)*monfn(t-p(3),force.force_fd);


r = cell(1,3);

r{1,1} = x{1};
r{1,2} = force;
r{1,3} = -p(2)*dforce;
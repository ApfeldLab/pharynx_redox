function r = delayd2fdp2(t,fd_cell,p,force)


force1 = eval_fd(t-p(3),force.force_fd,1);
dforce = exp(eval_fd(t-p(3),force.force_fd));
ddforce = force1.*dforce;

r = cell(1,3,3);
r(:) = {0};

r{1,1,3} = -dforce;
r{1,3,1} = -dforce;
r{1,3,3} = p(2)*ddforce;
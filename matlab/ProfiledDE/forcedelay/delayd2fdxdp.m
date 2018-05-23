function r = delayd2fdxdp(t,fd_cell,p,force)

r = cell(1,1,3);
r(:) = {0};

r{1,1,1} = 1;

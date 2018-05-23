function r = rossdfdx(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossdfdx
%
% The derivative of the Rossler equations with respect to x 
% at times t with parameters p. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



r = cell(3,3);
r(1:3,1:3) = {0};

r{1,2} = -1;
r{1,3} = -1;
r{2,1} = 1;
r{2,2} = p(1);
r{3,1} = eval_fd(t,fd_cell{3});
r{3,3} = -p(3) + eval_fd(t,fd_cell{1});


end
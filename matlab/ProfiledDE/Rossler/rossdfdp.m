function r = rossdfdp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossdfdp
%
% The derivative of the Rossler equations with respect to p 
% at times t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = cell2mat(eval_fdcell(t,fd_cell,0)');
% r = zeros(length(t),3,3);
% 
% %r(:,2,1) = y(:,2);
% r(:,2,1) = eval_fd(t,fd_cell{2});
% r(:,3,2) = 1;
% %r(:,3,3) = -y(:,3);
% r(:,3,3) = -eval_fd(t,fd_cell{3});

r = cell(3,3);
r(1:3,1:3) = {0};

r{2,1} = eval_fd(t,fd_cell{2});
r{3,2} = 1;
r{3,3} = -eval_fd(t,fd_cell{3});

end
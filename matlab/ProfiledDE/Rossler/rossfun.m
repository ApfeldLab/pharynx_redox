function r = rossfun(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossfun
%
% The Rossler equations in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y = cell2mat(eval_fdcell(t,fd_cell,0)');
% r = y;
% r(:,1) = -y(:,2) - y(:,3);
% r(:,2) = y(:,1) + p(1).*y(:,2);
% r(:,3) = p(2) + (y(:,1)-p(3)).*y(:,3);

y = eval_fdcell(t,fd_cell,0);
r = y;
r{1} = -y{2} - y{3};
r{2} = y{1} + p(1).*y{2};
r{3} = p(2) + (y{1}-p(3)).*y{3};

end
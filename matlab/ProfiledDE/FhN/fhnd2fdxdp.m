function r = fhnd2fdxdp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnd2fdxdp
%
% The second derivative of the FitzHugh-Nagumo equations with respect to 
% x and p at times t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = cell2mat(eval_fdcell(t,fd_cell,0)');
r = cell(2,2,3);
r(1:2,1:2,1:3) = {0};

%r(:,1,1,3) = 1-y(:,1).^2;
r{1,1,3} = 1 - eval_fd(t,fd_cell{1}).^2;
r{1,2,3} = 1;
r{2,1,3} = 1/p(3)^2;
r{2,2,2} = - 1/p(3);
r{2,2,3} = p(2)/p(3)^2;


end

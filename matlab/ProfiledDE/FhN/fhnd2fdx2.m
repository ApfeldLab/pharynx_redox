function r = fhnd2fdx2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnd2fdx2
%
% The second derivative of the FitzHugh Nagumo equations with respect 
% to the input x, at time t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = cell2mat(eval_fdcell(t,fd_cell,0)');
r = cell(2,2,2);
r(1:2,1:2,1:2) = {0};

%r(:,1,1,1) = -2*p(3)*y(:,1);
r{1,1,1} = -2*p(3)*eval_fd(t,fd_cell{1});
end
function r = fhndfdx(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhndfdx
%
% The derivative of the FitzHugh Nagumo equations with respect to x 
% at times t with parameters p. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = cell2mat(eval_fdcell(t,fd_cell,0)');
r = cell(2,2);
r(1:2,1:2) = {0};

%r(:,1,1) = (p(3)-p(3)*y(:,1).^2);
r{1,1} = p(3) - p(3)*eval_fd(t,fd_cell{1}).^2;
r{1,2} =  p(3);
r{2,1} = (-1/p(3));
r{2,2} = (-p(2)/p(3));

end
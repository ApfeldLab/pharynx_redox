function r = rossd3fdx3(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossd3fdy3
%
% The third derivative of the Rossler equations with respect to x.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = zeros(length(t),3,3,3,3);

r = cell(3,3,3,3);
r(1:3,1:3,1:3,1:3) = {0};

end
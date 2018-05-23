function r = rossd2fdx2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossd2fdx2
%
% The second derivative of the Rossler equations with respect 
% to the input x at time t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r = zeros(length(t),3,3,3);
%r(:,3,1,3) = 1;
%r(:,3,3,1) = 1;

r = cell(3,3,3);
r(1:3,1:3,1:3) = {0};
r{3,1,3} = 1;
r{3,3,1} = 1;

end
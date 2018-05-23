function r = rossd3fdx2dp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossd3fdx2dp
%
% The third derivative of the Rossler equations with respect to 
% x (twice) and p (once) at times t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = zeros(length(t),3,3,3,3);

r = cell(3,3,3,3);
r(1:3,1:3,1:3,1:3) = {0};

end
function r = rossd3fdxdp2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossd3fdxdp2
%
% The second derivative of the Rossler equations with respect to 
% x (once) and p (twice) at times t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = zeros(length(t),3,3,3,3);

r = cell(3,3,3,3);
r(1:3,1:3,1:3,1:3) = {0};

end
function r = rossd2fdxdp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rossd2fdxdp
%
% The second derivative of the Rossler equations with respect to 
% x and p at times t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = zeros(length(t),3,3,3);
%r(:,2,2,1) = 1;
%r(:,3,3,3) = -1;

r = cell(3,3,3);
r(1:3,1:3,1:3) = {0};

r{2,2,1} = 1;
r{3,3,3} = -1;

end
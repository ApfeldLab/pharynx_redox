function r = fhnd3fdxdp2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnd3fdxdp2
%
% The third derivative of the FitzHugh Nagumo equations with respect 
% to the input x and parameters p (twice) at time t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y = cell2mat(eval_fdcell(t,fd_cell,0)');
r = cell([2,2,3,3]);
r(1:2,1:2,1:3,1:3) = {0};

r{2,2,2,3} = 1/p(3)^2;
r{2,2,3,2} = 1/p(3)^2;
r{2,1,3,3} = - 2/p(3)^3;
r{2,2,3,3} = - 2*p(2)/p(3)^3;
end
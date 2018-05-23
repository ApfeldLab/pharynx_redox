function r = fhnd2fdp2(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnd2fdp2
%
% The second derivative of the FitzHugh Nagumo equations with respect 
% to the parameters p, at time t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);
r = cell(2,3,3);
r(1:2,1:3,1:3) = {0};
%r(:,1,1,1) = -2*p(3)*y(:,1);
r{2,3,1} =  -1/p(3)^2;
r{2,3,2} =  y{2}/p(3)^2;
r{2,1,3} =  -1/p(3)^2;
r{2,2,3} =  y{2}/p(3)^2;
r{2,3,3} =  - 2*(y{1}-p(1)+p(2)*y{2})/p(3)^3;
end
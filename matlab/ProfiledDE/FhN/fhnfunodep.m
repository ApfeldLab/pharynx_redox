function r = fhnfunodep(t,y,p,fd_obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnfunodep
%
% The FitzHugh-Nagumo equations in scalar form with perturbation given by
% fd_obj. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = y;
r(1) = p(3)*(y(1) - y(1).^3/3 + y(2)) + eval_fd(t,fd_obj);
r(2) = -(y(1) -p(1) + p(2)*y(2))/p(3);


end
function r = fhndfdp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhndfdp
%
% The derivative of the FitzHugh Nagumo equations with respect to p 
% at times t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = eval_fdcell(t,fd_cell,0);
r = cell(2,3);

r(1:2,1:3) = {0};

r{1,3} =  (y{1}-y{1}.^3/3+y{2});
r{2,1} = 1/p(3);
r{2,2} = (-y{2}/p(3));
r{2,3} = ((y{1}-p(1)+p(2)*y{2})/(p(3).^2));

end
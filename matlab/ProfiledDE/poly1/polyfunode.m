function r = polyfunode(t,y,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyfunode
%
% One dimensional polynomial function with forcing -- scalar form. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = 0*y;

for i = 1:(length(p)-1)
    r = r + p(i)*y.^i;
end 

r = r + p(length(p))*more.forcing(t,more.fs);

end
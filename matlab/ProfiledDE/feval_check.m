function [f] = feval_check(f,t,fd,p,extra)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% feval_check
%
% Evaluates a function f at times t with cell array fd and parameters p,
% may also have additional arguments extra. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(extra)
    f = feval(f,t,fd,p,extra);
else
    f = feval(f,t,fd,p);
end
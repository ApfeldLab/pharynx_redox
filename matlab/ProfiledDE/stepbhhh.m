%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,fs,iter,errcode]=stepbhhh(f,x0,f0,g0,d,MaxIters,varargin);
% STEPBHHH Computes an approximate minimum step length
%
% USAGE:
%  [s,fs,iters,errcode]=STEPBHHH(f,x0,f0,g0,d,MaxIters);
%
% INPUTS:
%        f : the objective function being minimized
%       x0 : the current value of the parameters
%       f0 : the value of f(x0)
%       g0 : the gradient vector of f at x0
%        d : the search direction       
% MaxIters : the maximum number of function evaluations allowed
%
% OUTPUTS:
%    s:  the optimal step in the direction d
%    fs: the value of f at x+s*d,
%    iter:  the number of iterations used
%    errcode: equals 1 if maximum iterations are exceeded
 
% STEPBHHH uses an algorithm based on one discussed in Berndt, et. al.,
% Annals of Economic and Social Measurement, 1974, pp. 653-665.
% This procedure specifies a cone of convergence in the plane
% defined by the direction vector, d, and the value of the objective
% function. The cone is defined by the lines through the origin
% (x,f(x)) with slopes (d'g)*delta and (d'g)*(1-delta).
% Delta must lie on (0,0.5).
% The procedure iterates until a point is found on the objective function
% that lies within the cone.
% In general, the wider the cone, the faster a "suitable" step size
% will be found.  If a trial point lies above the cone
% the step size will be increased and if it lies below the cone
% the step size is decreased. 
 
% ---------- INITIALIZATIONS ---------------------
  if nargin<6 | isempty(MaxIters), MaxIters=25; end
 % delta=optget('optstep','bhhhcone',0.0001);
  delta = 0.00000001;
  dg=-d'*g0;                   % directional derivative 
  tol1=dg*delta;
  tol0=dg*(1-delta);
  s=1;
  ds=1;
  errcode=0;
  
  % first bracket the cone
  for iter=1:MaxIters
    x=x0+s*d; fs=feval(f,x,varargin{:});
    temp=(f0-fs)/s;
    if temp<tol0
      ds=2*ds;
      s=s+ds;
    else break
    end
  end
  
  if tol0<=temp & temp<=tol1, return, end
  
  ds=ds/2;
  s=s-ds;
  it=iter+1;   
  
% then use bisection to get inside it
  for iter=it:MaxIters
    ds = ds/2;
    x=x0+s*d; fs=feval(f,x,varargin{:}); 
    temp=(f0-fs)/s;
    if     temp > tol1, s=s-ds; 
    elseif temp < tol0, s=s+ds;
    else return
    end  
  end
  errcode=1;
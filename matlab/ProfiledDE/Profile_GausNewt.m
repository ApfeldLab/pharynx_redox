function [newpars,DEfd,res0,Dres] = Profile_GausNewt(pars,lsopts,DEfd,fn,lambda,...
    Ycell,Tcell,wts,alg,lsopts2,fn_extras,pen,dpen,pen_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function Profile_GausNewt
%
% Gauss-Newton algorithm for use with ProfileErr_rep.
%
% INPUTS:
%
% pars  -   estimate of the parameters
%
% lsopts -  standard options for Optimization toolbox to be used in outer
%           Gaus-Newton scheme.
%
% parind -  an array of indices. The rows of parind provide the parameters
%           that apply to each replication of the experiment. 
%
% DEfd    - a cell array, each component contains a functional data 
%           object corresponding to the same component of the DE. This is
%           not entered, but assumed to be defined as a global variable.
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:       a handle to a function of the form f(t,y,p) giving the 
%               right hand side of the DE, should accept vector t, a cell 
%               array of function data objects, fd_cell, and output cell 
%               cell array of vectors, one for each component. 
%
%  fn.dfdx:     a handle to a function giving the derivative of fn.fn with 
%               respect to y. Output should be a cell array of vectors, the
%               first dimenion giving the components of fn.fn with
%               derivatives indexed by the second dimension. 
%
%  fn.dfdp:     a handle to a function giving the derivative of fn.fn with 
%               respect to p. Output should be a cell array of vectors, the
%               first dimenion giving the components of fn.fn with
%               derivatives indexed by the second dimension. 
%
%  fn.d2fdx2:   a handle to a function giving the Hessian of fn.fn with
%               respect to y. Output should be as above, with derivatives 
%               given in the second dimension of the cell array. 
% 
%  fn.d2fdxdp:  a handle to  a function giving the second derivative of 
%               fn.fn with respect to y and p. The first dimension sound
%               index components of fn.fn, with subsequent dimensions
%               indexing derivatives with respect to x and p in that order.
%
% lambda  - a vector of smoothing parameters corresponding to each
%           component of the DE.
%
% Ycell    - the Ycell of observed values of the DE, either matrix or cell
%           array.
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% lsopts2  - options for inner optimisation routine; may be empty.  
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% pen   - a penalty for the parameters; this should be a function that 
%         takes the parameters as arguments and additional input pen_extra
%
% dpen  - the derivative of the penalty with respect to parameters, this
%         should be a function that takes the parameters as arguments and
%         an additional input pen_extra
%
% pen_extras - additional input into penalty functions
%
% if wts is empty, the inverse standard deviation of each component
% is used as a weighting for each dimension. If lambda is a singleton, it 
% is multiplied by the average weight in each component. 
%
% OUTPUT
%
% newpars  - profiled estimates of the parameters
%
% DEfd     - estimated smooth of the differential equations at the profiled
%            parameters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           

if nargin<14
    pen_extras = [];
end

if nargin<13
    dpen = [];
end

if nargin<12
    pen = [];
end

if nargin<11
    fn_extras = [];
end

if nargin<10
    lsopts2 = [];
end

[res0,Dres,DEfd,dcdp] = ProfileErr(pars,DEfd,fn,lambda,Ycell,Tcell,wts,alg,...
    lsopts2,fn_extras,pen,dpen,pen_extras);

tres1 = res0;
tDres = Dres;
tDEfd = DEfd;
tdcdp = dcdp;

pars0 = pars;
pars1 = pars0;

F0 = sum(res0.^2);
F1 = F0;
iter = 0;
fundif = 1;
gradnorm0 = mean(res0'*Dres);
gradnorm1 = 1;
dcdp = [];

if isequal(lsopts.Display,'iter')
       header = sprintf(...
           '\n Iteration       steps    Residual   Improvement   Grad-norm     parameters');
       disp(header)
       formatstr = ' %5.0f       %5.0f   %13.6g %13.6g %12.3g     ';
%    disp(sprintf('\n'))
%    disp('iter    attempt        error         diff    grad norm     parameters')
end

while gradnorm1 > lsopts.TolFun && fundif > lsopts.TolFun
    
    
%    disp(pars0');
    iter = iter + 1;
    if iter > lsopts.MaxIter, break; end

%     if ~prod(+(eig(Dres'*Dres)<0)))
%         bit=min(eig(Dres'*Dres));
%         Dpars = (Dres+eye(size(Dres))*abs(2*bit))\res0;
%     else
        Dpars = Dres\res0;
%    end

    ntry = 0;
    while F1 >= F0 && norm(Dpars) > lsopts.TolFun && ntry < 10
    %    disp(F1)
        pars1 = pars0 - Dpars;
        [tres1,tDres,tDEfd,tdcdp] = ProfileErr(pars1,DEfd,fn,lambda,...
            Ycell,Tcell,wts,alg,lsopts2,fn_extras,pen,dpen,pen_extras,...
            dcdp,pars0);
        F1 = sum(tres1.^2);
        Dpars = Dpars/2;
        ntry = ntry + 1;
    end
    res1 = tres1;
    Dres = tDres;
    DEfd = tDEfd;
    dcdp = tdcdp;
    gradnorm1 = abs(mean(res1'*Dres));
    fundif = (F0-F1)/abs(F0);
    if isequal(lsopts.Display,'iter')
        iterstr = [sprintf(formatstr,iter,ntry,F1,fundif,gradnorm1),...
            num2str(pars1')];
        disp(iterstr);
        %        disp(num2str([iter, ntry, F1, fundif, gradnorm1,pars0']));
    end
    pars0 = pars1;
    res0 = res1;
    F0 = F1;
    gradnorm0 = gradnorm1;
end

newpars = pars0;

end
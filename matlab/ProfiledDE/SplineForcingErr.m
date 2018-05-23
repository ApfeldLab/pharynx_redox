function [f,J] = SplineForcingErr(bigcoef,basis_cell,basis_cellp,..
    whichindex,Ycell,Tcell,wts,lambda,lambdap,f_lfd,fn,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SplineForcingErr
%
% To be used with lsqnonlin; calculates forcing functions to make a
% differential equation agree with data
%
% INPUTS
%
% bigcoef  -  the concatenated co-efficients of the smooth and the forcing
%             function basis expanstion
%
% basis_cell  -  a cell-array of basis objects used for describing the
%             differential equation.
%
% basis_cellp -  a cell-array of basis objects used for describing the the
%             forcing functions. 
%
% whichindex - a vector giving the indices of the components to attach
%              forcing functions to. Should be the same dim as basis_cellp.
%
% Ycell   - observed values of the DE (matrix or cell array).
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
% 
% lambda      - a vector of smoothing co-efficients; will be transformed
%               into a vector if given as a scalar. 
%
% lambdap     - a vector of smoothing co-efficients for the forcing
%               functions.
%
% f_lfd       - an Lfd object specifying the smoothing penalty to penalize
%               the forcing functions. 
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:   a handle to a function of the form f(t,y,p) giving the right 
%           hand side of the DE, should accept vector t, cell array 
%           function data object, fd_cell, and output cell array of
%           vectors, one for each component. 
%
%  fn.dfdx: a handle to a function giving the derivative of fn.fn with 
%           respect to y. Output should be a cell array of vectors, the
%           first dimenion giving the components of fn.fn with derivatives
%           indexed by the second dimension. 
%
% pars  -   vector of inputs p into fn.fn. 
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
% 
% fn_extras   - object containing additional input to the right hand side 
%               of the differential equation. 
%
% if wts is empty, the inverse standard deviation of each component
% is used as a weighting. If lambda is a singleton, it is multiplied by
% the average weight in each dimension. 
%
% OUTPUT:
%
% f  -  vector residuals of the fit, and of the penalty at quadrature
%       points.
%
% J  -  derivatives of f with respect to coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<14
    fn_extras = [];
end

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

lambda = sqrt(lambda);

if length(lambdap)==1
    lambdap = lambdap*ones(size(whichindex));
    lambdap = reshape(lambdap,numel(lambdap)),1);
end

% if length(lambdap)<length(basis_cell))
%     lambdap2 = zeros(length(basis_cell),1);
%     lambdap2(whichindex) = lambdap;
%     lambdap = lambdap2;
% end

if length(whichindex) ~= length(basis_cellp)
    error('Bases for forcing functions do not agree with whichindex.')
end

wi = zeros(size(Ycell,2));
wi(whichindex) = 1;


% Calculate coefsizes

coefsizes = 0;
for i = 1:numel(basis_cell)
   coefsizes = coefsizes + getnbasis(basis_cell{i});
end

% Re-pass parameter-coefficient matrix:

pcoef = bigcoef(1:coefsizes);
newcoefs = bigcoef((coefsizes+1):length(bigcoef));


% Create the appropriate functional data object

DEfd_cell = Make_fdcell(newcoefs,basis_cell);
pDEfd_cell = Make_fdcell(pcoef,basis_cellp);

% Matrix of errors

predns = reshape(eval_fdcell(Tcell,DEfd_cell,0),size(Ycell,1),size(Ycell,2));

% Matrix of penalty at quadrature points

quadvals = getquadvals(basis_cell{1});
qvals = sqrt(quadvals(:,2));
qpredns = cell2mat(eval_fdcell(quadvals(:,1),DEfd_cell,0));
dpredns = cell2mat(eval_fdcell(quadvals(:,1),DEfd_cell,1));

ppredns = zeros(size(quadvals,1),size(Ycell,2));
ppredns(:,whichindex) = eval_fdcell(quadvals(:,1),pDEfd_cell,0);

fvals = cell2mat(feval_check(fn.fn,quadvals(:,1),qpredns,pars,fn_extras));
dfvals = feval_check(fn.dfdx,quadvals(:,1),qpredns,pars,fn_extras);

pen = reshape(dpredns - fvals - ppredns);

pen = kron(qvals,lambda).*pen;

pen = reshape(pen,size(pen,1)*size(pen,2),1);



% Standard penalties for forcing function


pen2 = eval_fdcell(quadvals(:,1),pDEfd_cell,f_lfd);
pen2 = kron(qvals,lambdap).*pen2;
pen2 = reshape(pen2,size(pen2,1)*size(pen2,2),1);

dpen = eval_basis_cell(quadvals(:,1),pDEfd_cell,pen_lfd);

for i = 1:length(pDEfd_cell)
    dpen{i} = dpen{i}.*kron(qvals,lambdap);
end

dpen = mattdiag_cell(dpen,0);


% Now to calculate a Jacobian

if nargout > 1

    [err,Zmat] = djdc_pts(fd_cell,Ycell,Tcell,wts);
    
    Dphimat_cell = getvalues_cell(basis_cell,1);
    phimat_cell = getvalues_cell(basis_cell,0);

    pphimat_cell = getvalues_cell(basis_cellp,0);
    
% Derviatives of penalty wrt coefs:

    
    dpen2 = cell(length(phimat_cell));
    
    for i = 1:size(Ycell,2)
        Dphimat_cell{i} = spdiags(qvals,0,m,m)*...
            sparse(lambda(i)*Dphimat_cell{i});
        for j = 1:size(Ycell,2)
              dpen2{i,j}=spdiags(dfvals{i,j}.*qvals,0,m,m)*...
                sparse(lambda(i)*phimat_cell{j});
        end
        
    end

    dpen1 = mattdiag_cell(Dphimat_cell,0);
    dpen2 = cell2mat(dpen2);

    J2 = [];
    wdiag=diag(wi);
    k = 0;
    for i = 1:size(Ycell,2)
       if wdiag(i)
           k = k+1;
           J2 = [J2; kron(wi(:,i),lambda(i)*diag(qvals)*pphimat_cell{k})];
       end
    end
    
    
    J1 = [-Zmat; (dpen1-dpen2)];
    J2 = [zeros(size(Zmat);J2];
  
    J = -[J2 -J1];    
    
    J = [J; [lambdap*dpen zeros(1,size(J,2)-length(dpen))]];
else 
    
    err = djdc_pts(fd_cell,Ycell,Tcell,wts);
    
end

    % Putting the error vector together. 
    
    f = [err; pen; pen2];

end

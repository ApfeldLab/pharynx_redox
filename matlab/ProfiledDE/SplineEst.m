function [DEfd_cell,resnorm] = SplineEst(fn,Tcell,Ycell,pars,...
    knots_cell,wts,lambda,lambda_first,rough_ord,alg,lsopts,...
    DEfd_cell,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SplineEst
%
% A function to create a spline smooth of a Differential Equation. 
%
% INPUTS:
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
% Tcell - vector of times at which the DE has been evaluated. May be either
%         a vector or cell array. 
%
% Ycell  - observations of the DE at times given in Tcell. May be either a
%         matrix (if Tcell a vector) or a cell array. 
%
% pars  - parameters p for the functions fn.fn and fn.dfdx.
%
% knots_cell - a cell array giving the knots to use in a basis expansion.
%              each component of the array corresponds to a dimension of
%              the DE. 
%
% wts        - a vector of weights to give the squared error of each
%              component in the fit. May be empty
%
% lambda     - a vector of smoothing parameters.
%
% lambda_first - smoothing parameters for a roughness penalty. 
%
% rough_ord - order of roughness penalty. 
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% lsopts     - options for the least squares routine. May be empty. 
%
% DEfd_cell  - an existing approximate smooth is present. In this case no
%              attempt is made to smooth the data with a roughness penalty.
%
% fn_extras   - object containing additional input to the right hand side 
%               of the differential equation. 
%
% if wts is empty, the fit is weighted by the inverse of the standard
% of each component of the Ycell. If lambda is a singleton it is turned 
% into a vector by multiplying it by wts. 
%
% OUTPUT:
%
% DEfd_cell  - a cell array. Each component contains a functional data
%              object that is the fit for that component of the DE. 
%
% resnorm    - the residual norm of the fit. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<13
    fn_extras = [];
end

if nargin<12
    DEfd_cell = [];
end

% Turn Ycell into a cell array if needed

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

if isempty(alg) 
    alg = ones(length(Ycell),1);
end

% Set up quadrature points;
if isempty(DEfd_cell)
    nquad = 21;

    if ~iscell(knots_cell)
        knots = cell(1,length(Ycell));
        knots(:) = {knots_cell};
        knots_cell = knots;
    end
    
    if size(knots_cell,1)>1
        knots_cell = knots_cell';
    end
    Tcell2 = cell2mat(knots_cell);
    
    range = [min(Tcell2) max(Tcell2)];

    basis_cell = cell(1,length(Ycell));

    nbasis = zeros(length(Ycell),1);
    norder = 6;


    bigknots = knots_cell{1};
    nbasis(1) = length(knots_cell{1}) + norder - 2;

    for i = 2:length(Ycell)
        bigknots = [bigknots knots_cell{i}];
        nbasis(i) = length(knots_cell{i}) + norder -2;
    end

    quadvals = MakeQuadPoints(bigknots,nquad);


    Lfd_cell = cell(1,length(Ycell));

    for i = 1:length(Ycell)
        basis_cell{i} = MakeBasis(range,nbasis(i),norder,knots_cell{i},...
            quadvals,1);
        Lfd_cell{i} = fdPar(basis_cell{i},rough_ord,lambda_first);
    end

    %DEfd_cell = data2fd_cell(basis_cell,Ycell,Tcell);

    DEfd_cell = smoothfd_cell(Ycell,Tcell,Lfd_cell)';
else
    basis_cell = getcellbasis(DEfd_cell);
end

coefs = getcellcoefs(DEfd_cell);

[newcoefs,resnorm] = lsqnonlin(@SplineCoefErr,coefs,[],[],lsopts,...
    basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,fn_extras);

DEfd_cell = Make_fdcell(newcoefs,basis_cell);


% figure(1)
% for i = 1:length(Ycell))
%     subplot(length(Ycell),1,i);
%     plotfit_fd(Ycell{i},Tcell{i},DEfd_cell{i});
% end


end

function d2G = make_d2gdc2(fd_cell,fn,Tcell,wts,lambda,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdc2
%
% The second derivative of a spline objective function with respect to 
% the coefficients.
%
% INPUTS
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:     a handle to a function of the form f(t,y,p) giving the right 
%             hand side of the DE, should accept vector t, cell array 
%             function data object, fd_cell, and output cell array of
%             vectors, one for each component. 
%
%  fn.dfdx:   a handle to a function giving the derivative of fn.fn with 
%             respect to y. Output should be a cell array of vectors, the
%             first dimenion giving the components of fn.fn with 
%             derivatives indexed by the second dimension. 
%
%  fn.d2fdx2: a handle to a function giving the Hessian of fn.fn with 
%             respect to y. Output should be as above, with derivatives 
%             given in the second and third dimension of the cell array. 
% 
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. If
%           empty, all weights are assumed equal.
%
% lambda  - a vector of smoothing parameters corresponding to each
%           component of the DE.
%
% pars    - vector of the current parameters in the penalty.
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% OUTPUT:
%
% d2G   -   the second derivative matrix of the spline criterion with
%           respect to coefficients. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<8
    fn_extras = [];
end

% Convert wts, lambda and Ycell to the right format:

Ycell = eval_fdcell(Tcell,fd_cell,0);
[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(Ycell),1);
end

% disp('setting up')

basis_cell = getcellbasis(fd_cell);

quads = getquadvals(basis_cell{1});

n = size(quads,1);

quadpts = quads(:,1);
quadvals = quads(:,2);

dpath = eval_fdcell(quadpts,fd_cell,1);

fpath = feval_check(fn.fn,quadpts,fd_cell,pars,fn_extras);
dfpath = feval_check(fn.dfdx,quadpts,fd_cell,pars,fn_extras);
d2fpath = feval_check(fn.d2fdx2,quadpts,fd_cell,pars,fn_extras);


% disp('calculating penalty weight values')


wt2 = cell(length(basis_cell),length(basis_cell));
wt2(1:length(basis_cell),1:length(basis_cell)) = {0};

for i = 1:length(basis_cell)
    for j = 1:length(basis_cell)
        for k = 1:length(basis_cell)
            wt2{i,j} = wt2{i,j} + lambda(k)*(d2fpath{k,i,j}.*(fpath{k}-dpath{k}) + ...
                dfpath{k,i}.*dfpath{k,j});
        end
    end
end


% disp('putting it together')

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);

Zmat_cell = eval_basis_cell(Tcell,basis_cell,0);

H = cell(length(phimat_cell));

for i = 1:length(basis_cell)
    for j = 1:length(basis_cell)
        wt2{i,j} = wt2{i,j}.*quadvals;
        H{i,j} = phimat_cell{i}'*spdiags(wt2{i,j},0,n,n)*phimat_cell{j} - ... 
            lambda(i)*Dphimat_cell{i}'*spdiags(dfpath{i,j}.*quadvals,0,n,n)*phimat_cell{j} - ...
            lambda(j)*phimat_cell{i}'*spdiags(dfpath{j,i}.*quadvals,0,n,n)*Dphimat_cell{j};
    end

    H{i,i} = H{i,i} + ...
        lambda(i)*Dphimat_cell{i}'*spdiags(quadvals,0,n,n)*Dphimat_cell{i};
    
    if ~isempty(Tcell{i})
        m = length(wts{i});
        H{i,i} = H{i,i} + Zmat_cell{i}'*spdiags(wts{i},0,m,m)*Zmat_cell{i};
    end
end        

% disp('making the final matrix')

d2G = sparse(cell2mat(H));

end
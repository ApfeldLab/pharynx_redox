function d3G = make_d3gdc2dp(fd_cell,fn,lambda,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d3jdc2dp
%
% The third derivative of a spline objective function with respect to 
% the coefficients (twice) and the parameters (once).
%
% INPUTS
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
%  fn.fn:       a handle to a function of the form f(t,y,p) giving the
%               right hand side of the DE, should accept vector t, cell
%               array of function data object, fd_cell, and output cell 
%               array vectors, one for each component. 
%
%  fn.dfdx:     a handle to a function giving the derivative of fn.fn with 
%               respect to y. Output should be a cell array of vectors, the
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
%  fn.d3fdx2dp: a handle to a function giving the third derivative of fn.fn
%               with respect to y twice (second and third dimensions) and p
%               (fourth dimension). 
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
% d3G   -   An cell array of the same length as the parameters. Each element
%           contains the third derivative matrix of the spline criterion
%           with respect to the coefficients and one parameter. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6 
    fn_extras = [];
end

% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(fd_cell),1);
end

% disp('setting up')

basis_cell = getcellbasis(fd_cell);

quads = getquadvals(basis_cell{1});

n = size(quads,1);

quadpts = quads(:,1);
quadvals = quads(:,2);

dYcell = eval_fdcell(quadpts,fd_cell,1);

fYcell = feval_check(fn.fn,quadpts,fd_cell,pars,fn_extras);
fn.dfdxYcell = feval_check(fn.dfdx,quadpts,fd_cell,pars,fn_extras);
fn.dfdpYcell = feval_check(fn.dfdp,quadpts,fd_cell,pars,fn_extras);
fn.d2fdx2Ycell = feval_check(fn.d2fdx2,quadpts,fd_cell,pars,fn_extras);
fn.d2fdxdpYcell = feval_check(fn.d2fdxdp,quadpts,fd_cell,pars,fn_extras);
d3fYcell = feval_check(fn.d3fdx2dp,quadpts,fd_cell,pars,fn_extras);


% disp('calculating penalty weight values')

wt2 = cell(length(basis_cell),length(basis_cell),length(pars));
wt2(1:length(basis_cell),1:length(basis_cell),1:length(pars)) = {0};

for i = 1:length(basis_cell)
    for j = 1:length(basis_cell)
        for k = 1:length(basis_cell)
            for l = 1:length(pars) 
                wt2{i,j,l} = wt2{i,j,l} + lambda(k)*(...
                    d3fYcell{k,i,j,l}.*(fYcell{k}-dYcell{k}) + ...
                    fn.d2fdx2Ycell{k,i,j}.*fn.dfdpYcell{k,l} + ...
                    fn.d2fdxdpYcell{k,i,l}.*fn.dfdxYcell{k,j} + ...
                    fn.d2fdxdpYcell{k,j,l}.*fn.dfdxYcell{k,i});
            end
        end
    end
end

% disp('putting it together')

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);

H = cell([length(phimat_cell) length(phimat_cell) length(pars)]);
d3G = cell(length(pars),1);


for l = 1:length(pars)
    for i = 1:length(basis_cell)
        for j = 1:length(basis_cell)
            wt2{i,j,l} = wt2{i,j,l}.*quadvals;
            H{i,j,l} = phimat_cell{i}'*spdiags(wt2{i,j,l},0,n,n)*phimat_cell{j} - ...
                lambda(i)*Dphimat_cell{i}'*spdiags(fn.d2fdxdpYcell{i,j,l}.*quadvals,0,n,n)*phimat_cell{j} - ...
                lambda(j)*phimat_cell{i}'*spdiags(fn.d2fdxdpYcell{j,i,l}.*quadvals,0,n,n)*Dphimat_cell{j};
        end
    end

    d3G{l} = sparse(cell2mat(H(:,:,l)));
end

end
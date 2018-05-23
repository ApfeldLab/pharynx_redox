function d3G = make_d3gdc3(fd_cell,fn,lambda,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d3jdc3
%
% The third derivative of a spline objective function with respect to 
% the coefficients.
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
%  fn.d3fdx3:   a handle to a function giving the third derivative of fn.fn
%               2nd dimnesion) thrice with respect to y. 
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
% d3G   -   An cell array of the same length as the co-efficients. Each element
%           contains the third derivative matrix of the spline criterion
%           with respect to the coefficients. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    fn_extras = [];
end

% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(path),1);
end

% disp('setting up')

basis_cell = getcellbasis(fd_cell);

quads = getquadvals(basis_cell{1});

n = size(quads,1);

quadpts = quads(:,1);
quadvals = quads(:,2);

dpath = eval_fdcell(quadpts,fd_cell,1);

fpath = feval_check(fn.fn,quadpts,fd_cell,pars,fn_extras);
dfdxpath = feval_check(fn.dfdx,quadpts,fd_cell,pars,fn_extras);
d2fdx2path = feval_check(fn.d2fdx2,quadpts,fd_cell,pars,fn_extras);
d3fdx3path = feval_check(fn.d3fdx3,quadpts,fd_cell,pars,fn_extras);

% disp('calculating penalty weight values')

wt2 = cell(length(basis_cell),length(basis_cell),length(basis_cell));
wt2(1:length(basis_cell),1:length(basis_cell),1:length(basis_cell)) = {0};

for i = 1:length(basis_cell)
    for j = 1:length(basis_cell)
        for k = 1:length(basis_cell)
            for l = 1:length(basis_cell)
                wt2{i,j,l} = wt2{i,j,l} + lambda(k)*(...
                    d3fdx3path{k,i,j,l}.*(fpath{k}-dpath{k}) + ...
                    d2fdx2path{k,i,j}.*dfdxpath{k,l} + ...
                    d2fdx2path{k,i,l}.*dfdxpath{k,j} + ...
                    d2fdx2path{k,j,l}.*dfdxpath{k,i});
            end
        end
    end
end

% disp('putting it together')

m = 0;
for(i = 1:length(basis_cell))
    m = m + getnbasis(basis_cell{i});
end

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);

H = cell([length(basis_cell) length(basis_cell)]);
d3G = cell(m,1);


here = 0;
for l = 1:length(basis_cell)
    nbas = getnbasis(basis_cell{l});
    for k = 1:nbas
        s = here+k;
        for i = 1:length(basis_cell)
            for j = 1:length(basis_cell)
                wt3 = wt2{i,j,l}.*quadvals.*phimat_cell{l}(:,k);
                H{i,j} = phimat_cell{i}'*spdiags(wt3,0,n,n)*phimat_cell{j} - ...
                    lambda(i)*Dphimat_cell{i}'*...
                        spdiags(d2fdx2path{i,j,l}.*quadvals.*...
                        phimat_cell{l}(:,k),0,n,n)*phimat_cell{j} - ...
                    lambda(j)*phimat_cell{i}'*...
                        spdiags(d2fdx2path{j,i,l}.*quadvals.*...
                        phimat_cell{l}(:,k),0,n,n)*Dphimat_cell{j} - ...
                    lambda(l)*phimat_cell{i}'*...
                        spdiags(d2fdx2path{l,j,i}.*quadvals.*...
                        Dphimat_cell{l}(:,k),0,n,n)*phimat_cell{j};
            end
        end
        d3G{s} = sparse(cell2mat(H));
    end
    here = here+nbas;
end

    
end
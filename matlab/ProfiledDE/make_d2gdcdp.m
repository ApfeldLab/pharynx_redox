function d2G = make_d2gdcdp(fd_cell,fn,lambda,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdcdp
%
% The second derivative of a spline objective function with respect to 
% the coefficients and parameters in a penalty function.
%
% INPUTS
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:      a handle to a function of the form f(t,y,p) giving the right 
%              hand side of the DE, should accept vector t, cell array 
%              function data object, fd_cell, and output cell array of
%              vectors, one for each component. 
%
%  fn.dfdx:    a handle to a function giving the derivative of fn.fn with 
%              respect to y. Output should be a cell array of vectors, the
%              first dimenion gives the components of fn.fn with 
%              derivatives indexed by the second dimension. 
%
%  fn.dfdp:    a handle to a function giving the derivative of fn.fn with 
%              respect to p. Output should be a cell array of vectors, the
%              first dimenion gives the components of fn.fn with 
%              derivatives indexed by the second dimension. 
% 
%  fn.d2fdxdp: a handle to a function giving the second derivative of fn.fn 
%              with respect to y and p. Should return a cell array of
%              vectors with dimensions indexing components of the system,
%              derivatives with respect to y and derivatives with respect
%              to p in that order. 
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
%           respect to coefficients and parameters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    fn_extras = [];
end

% disp('setting up')

if length(lambda)==1
    lambda = lambda*ones(size(fd_cell));
elseif  ~prod(+((size(lambda)==size(fd_cell))))
    disp(size(lambda))
    disp(size(fd_cell))
    error('Wrong number of lambdas');
end

if isempty(alg) 
    alg = ones(length(fd_cell),1);
end


basis_cell = getcellbasis(fd_cell);

quads = getquadvals(basis_cell{1});

n = size(quads,1);

quadpts = quads(:,1);
quadvals = quads(:,2);

dpath = eval_fdcell(quadpts,fd_cell,1);

fpath = feval_check(fn.fn,quadpts,fd_cell,pars,fn_extras);
dfpath = feval_check(fn.dfdx,quadpts,fd_cell,pars,fn_extras);
dfpathp = feval_check(fn.dfdp,quadpts,fd_cell,pars,fn_extras);
d2fpath = feval_check(fn.d2fdxdp,quadpts,fd_cell,pars,fn_extras); 


%disp('weight values')


wt2 = cell(length(basis_cell),length(pars));
wt2(1:length(basis_cell),1:length(pars)) = {0};

for j = 1:length(basis_cell)
    for i = 1:length(basis_cell)
        for k = 1:length(pars)
            wt2{i,k} = wt2{i,k} + lambda(j)*(d2fpath{j,i,k}.*(fpath{j}-dpath{j}) + ...
                dfpath{j,i}.*dfpathp{j,k});
        end
    end
end


% disp('putting it together')

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);

H = cell(length(phimat_cell),length(pars));


for i = 1:length(basis_cell)
    for j = 1:length(pars)
        H{i,j} = phimat_cell{i}'*(quadvals.*wt2{i,j}) - ...
            lambda(i)*Dphimat_cell{i}'*(quadvals.*dfpathp{i,j});
end

% disp('Making the matrix')

d2G = cell2mat(H);

end
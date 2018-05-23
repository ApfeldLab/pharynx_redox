function d3G = make_d3gdcdp2(fd_cell,fn,lambda,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d3jdcdp2
%
% The third derivative of a spline objective function with respect to 
% the coefficients (once) and parameters (twice) in a penalty function.
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
%              first dimenion giving the components of fn.fn with
%              derivatives indexed by the second dimension. 
%
%  fn.dfdp:    a handle to a function giving the derivative of fn.fn with 
%              respect to p. Output should be a cell array of vectors, the
%              first dimenion giving the components of fn.fn with
%              derivatives indexed by the second dimension.
% 
%  fn.d2fdxdp: a handle to  a function giving the second derivative of
%              fn.fn with respect to y (second dimension) and p (third 
%              dimension).  
%
%  fn.d2fdp2:   a handle to a function giving the second derivative of
%               fn.fn with respect to p (second and third dimensions). 
%
%  fn.d3fdxdp2: a handle to a function giving the third derivative of fn.fn
%               with respect to y (second dimension) and p twice (third and 
%               fourth dimensions). 
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
% d3G   -   a cell array of the same length as the parameters each
%           component containing a third derivative defined with respect to
%           the co-efficients, the parameters and one of the parameters
%           taken again. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    fn_extras = [];
end

% disp('setting up')

if length(lambda)==1
    lambda = lambda*ones(size(fd_cell));
elseif ~(size(lambda)==size(fd_cell))
    error('Wrong number of lambdas');
end

if isempty(alg) 
    alg = ones(length(path),1);
end


basis_cell = getcellbasis(fd_cell);

quads = getquadvals(basis_cell{1});

n = size(quads,1);

quadpts = quads(:,1);
quadvals = quads(:,2);

dpath = eval_fdcell(quadpts,fd_cell,1);

fpath = feval_check(fn.fn,quadpts,fd_cell,pars,fn_extras);
dfdxpath = feval_check(fn.dfdx,quadpts,fd_cell,pars,fn_extras);
dfdppath = feval_check(fn.dfdp,quadpts,fd_cell,pars,fn_extras);
d2fdxdppath = feval_check(fn.d2fdxdp,quadpts,fd_cell,pars,fn_extras); 
d2fdp2path = feval_check(fn.d2fdp2,quadpts,fd_cell,pars,fn_extras);
d3fpath = feval_check(fn.d3fdxdp2,quadpts,fd_cell,pars,fn_extras);

%disp('weight values')

wt2 = cell(length(basis_cell),length(pars),length(pars));
wt2(1:length(basis_cell),1:length(pars),1:length(pars)) = {0};

for k = 1:length(basis_cell)
    for i = 1:length(basis_cell)
        for j = 1:length(pars)
            for l = 1:length(pars)
                wt2{i,j,l} = wt2{i,j,l} + lambda(k)*(d3fpath{k,i,j,l}.*(fpath{k}-dpath{k}) + ...
                    d2fdxdppath{k,i,j}.*dfdppath{k,l} + ...
                    d2fdxdppath{k,i,l}.*dfdppath{k,j} + ...
                    dfdxpath{k,i}.*d2fdp2path{k,j,l});
            end
        end
    end
end


% disp('putting it together')

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);

H = cell(length(phimat_cell),length(pars));

d3G = cell(length(pars),1);

for l = 1:length(pars)
    for i = 1:length(basis_cell)
        for j = 1:length(pars)
            H{i,j} = (phimat_cell{i}'*(quadvals.*wt2{i,j,l}) - ...
                lambda(i)*Dphimat_cell{i}'*(quadvals.*d2fdp2path{i,j,l}));
        end

        % disp('Making the matrix')
        d3G{l} = sparse(cell2mat(H));
    end

end


end
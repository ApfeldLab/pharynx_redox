function V = make_varC(fd_cell,fn,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdc2
%
% Approximate covariance for coefficients in a spline approximation due to
% stochastic re-sampling. 
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
% lambda  - a vector of estimated stochastic variancesSEI_fd
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


if nargin<5
    fn_extras = [];
end

% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(fd_cell),1);
end

% disp('setting up')

basis_cell = getcellbasis(fd_cell);

quads = getquadvals(basis_cell{1});
coefs = getcellcoefs(fd_cell);

n = size(quads,1);

quadpts = quads(:,1);
quadvals = quads(:,2);


dfpath = feval_check(fn.dfdx,quadpts,fd_cell,pars,fn_extras);
lambda = feval_check(fn.procvar,quadpts,fd_cell,pars,fn_extras);


% disp('calculating penalty weight values')


wt2 = cell(length(basis_cell),length(basis_cell));
wt2(1:length(basis_cell),1:length(basis_cell)) = {0};
wt1 = wt2;


for i = 1:length(basis_cell)
    for j = 1:length(basis_cell)
        for k = 1:length(basis_cell)
            for l = 1:length(basis_cell)
                wt2{i,j} = wt2{i,j} + lambda{k,l}.*dfpath{k,i}.*dfpath{l,j};
            end
            wt1{i,j} = wt1{i,j} + lambda{i,k}.*dfpath{k,j};
        end
    end
end


% disp('putting it together')

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);

H = cell(length(phimat_cell));



for i = 1:length(basis_cell)
    for j = 1:length(basis_cell)
        wt2{i,j} = wt2{i,j}.*quadvals;
        H{i,j} = phimat_cell{i}'*spdiags(wt2{i,j},0,n,n)*phimat_cell{j} - ... 
            Dphimat_cell{i}'*spdiags(wt1{i,j}.*quadvals,0,n,n)*phimat_cell{j} - ...
            phimat_cell{i}'*spdiags(wt1{j,i}.*quadvals,0,n,n)*Dphimat_cell{j};
        H{i,j} = H{i,j} + Dphimat_cell{i}'*spdiags(lambda{i,j}.*quadvals,0,n,n)*Dphimat_cell{i};
    end
end        

% disp('making the final matrix')

d2G = sparse(cell2mat(H));

V = inv(0.99*d2G + 0.01*diag(diag(d2G)));

end
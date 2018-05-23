function [coefs,fd_cell,dcdy,df] = genlin_smooth(Ycell,Tcell,wts,basis_cell,...
    lambda,pars,alg,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genlin_smooth
%
% produces an exactly calculated smooth with penalty corresponding to a 
% forced, linear differential equation. 
%
% INPUTS
%
% Ycell   - observed values of the DE (matrix or cell array).
% 
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. If a
%           link function is used, Tcell must be a vector and different
%           measurements should be selected by setting some of the entries
%           in wts to zero. 
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% basis_cell - a cell-array, each component should contain a basis for
%                the corresponding dimension of the DE. Assumes quadrature
%                values have been attached to each. 
% 
% lambda - vector of smoothing parameters for each dimension
%
% pars  -  vector of inputs parameters
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% more  -   a struct specifying the entries in the linear differential
%           equation \dot(x) = Ax + Bu. It may have elements
%           
%           - mat:  the default matrix for the linear term A. This
%                   defaults to all zero.
%           - sub:  a k-by-2 array specifying the position in more.mat
%                   that correspond to parameters. By default this is all
%                   of them taken row-wise. 
%           - force: a cell array of forcing functions; either functions
%                    that take arguments t (plus an optional extra argument) 
%                    and return values at t, or functional data objects,
%                    assumed to be empty if not given
%           - force_mat: default matrix for B, again assumed to be all
%                        zeros
%           - force_sub: an l-by-2 mapping from unused values of pars to
%                        entries of B. Assumed to cover all of B by rows if
%                        not provided
%
% OUTPUT:
%
% coefs - a set of co-efficients describing the smooth
%
% fd_cell - a cell array of functional data objects describing the smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<8, more = []; end
if nargin<7, alg = []; end

% Convert wts, lambda and Ycell to the right format:

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

% By default, system has no algebraic components. 

if isempty(alg), alg = ones(length(Ycell),1); end

% Now put together the linear matrix and forcing terms

n = size(basis_cell,2);
more = checkmore(more,n);
quadvals = getquadvals(basis_cell{1});

pmat = more.mat;
pmat(sub2ind([n n],more.sub(:,1),more.sub(:,2))) = pars(1:size(more.sub,1));

if ~isempty(more.force)
    t = quadvals(:,1);

    fs = zeros(length(t),length(more.force));
    for i = 1:length(more.force)
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force_input);
        end
    end

    b = more.force_mat;
    b(sub2ind(size(b),more.force_sub(:,1),more.force_sub(:,2))) = ...
        pars((size(more.sub,1)+1):(size(more.sub,1)+size(more.force_sub,1)));

    b = reshape(fs*b',size(fs,1)*size(b,1),1);
else 
    b = zeros(n*size(quadvals,1),1);
end

% Put basis expressions together

y = cell2mat(Ycell');

Dphimat_cell = getvalues_cell(basis_cell,alg);
phimat_cell = getvalues_cell(basis_cell,0);


Zmat_cell = eval_basis_cell(Tcell,basis_cell,0);

Zmat = mattdiag_cell(Zmat_cell,0);

phimat_cell = repmat(phimat_cell,n,1);

for i = 1:n
    for j = 1:n
        phimat_cell{i,j} = pmat(i,j)*phimat_cell{i,j};
    end
end

Dphimat = mattdiag_cell(Dphimat_cell,0);
phimat = cell2mat(phimat_cell);

qvals = kron(quadvals(:,2),lambda');
wts = cell2mat(wts');

% Evaluate a large linear system

coefs = (Zmat'*diag(wts)*Zmat + (Dphimat-phimat)'*diag(qvals)*(Dphimat-phimat))\...
    (Zmat'*diag(wts)*y + (Dphimat-phimat)'*diag(qvals)*b);

if nargout>1
    fd_cell = Make_fdcell(coefs,basis_cell);

    if nargout > 2
        dcdy =  (Zmat'*diag(wts)*Zmat + (Dphimat-phimat)'*diag(qvals)*(Dphimat-phimat))\...
            (Zmat'*diag(wts));

        if nargout > 3
            df = trace((eye(size(Zmat,1)) - Zmat*dcdy)'*(eye(size(Zmat,1))-Zmat*dcdy));
        end

    end

end
function mats = eval_basis_cell(t,basis_cell,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function eval_basis_cell
%
% creates a cell array by evaluating a cell array of basis cells at times t
% (either a cell array or a fixed vector of common time points for all 
% elements of basis_cell) with derivative degree k.
%
% INPUTS:
%
% t - a vector of time points to evaluate the bases. May also be a cell
%     array of vectors, one element corresponding to one basis. 
%
% basis_cell - a cell array of basis objects
%
% k - the derivative of the basis functions to be evaluated. If a
%     the same derivative is evaluated accross all bases. If an array, the
%     each basis is evaluated at the derivative given in the corresponding
%     entry of k
%
% OUTPUT
%
% mats - a cell array of matrices corresponding to the evaluation of the
%        basis objects. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mats = cell(size(basis_cell));

t2 = t;

if length(k)==1  
    k = k*ones(size(basis_cell));
end
if prod(+(size(k)==size(basis_cell)))

    for i = 1:numel(basis_cell)
        if iscell(t)
            t2 = t{i};
        end
        if ~isempty(t2)
            mats{i} = eval_basis(t2,basis_cell{i},k(i));
        else
            mats{i} = zeros(0,getnbasis(basis_cell{i}));
        end
    end
else
    error('k is of the wrong dimension')
end
function [wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weightslambdasYcell
%
% This function converts objects for wts, lambda and Ycell into cell arrays
% if they are given in other forms. 
%
% INPUTS
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. If
%           empty, it is assumed all ones. 
%
% lambda - vector of smoothing parameters for each dimension, if a
%          singleton, lambda is multiplied by the average weight in
%          each dimension. 
%
% Ycell - a cell array giving the observed values of the Ycell in each
%        dimension
%
% OUTPUTS
%
% wts   - a cell array of weights for each observations
%
% lamda - a vector of smoothing parameters
%
% Ycell  - a cell array giving the observations in each dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Convert Ycell if necessary

if ~iscell(Ycell)
    Ycell_cell = cell(1,size(Ycell,2));
    for i = 1:size(Ycell,2)
        Ycell_cell{i}= Ycell(:,i);
    end
    Ycell = Ycell_cell;
end


% Now deal with a singleton weight value - transform into a vector

if ~iscell(wts)
    if isempty(wts)
        wts = ones(size(Ycell));
    elseif length(wts)==1
        wts = wts*ones(size(Ycell));
    elseif size(wts,1)==1
        wts = repmat(wts,size(Ycell,1),1);
    end

    if size(wts)~=size(Ycell)
        error('Weights must have same dimension as Ycell.');
    end

    wts2 = cell(size(Ycell));

    for i = 1:numel(Ycell)
        wts2{i} = wts(i)*ones(size(Ycell{i}));
    end
    
    wts = wts2;

else
    for i=1:numel(Ycell)
        if length(wts{i})~=length(Ycell{i})
            error('Weights must have same dimension as Ycell.');
        end
    end
end
% Now turn lambda into a vector if it was a singleton

if length(lambda) == 1
    lambda = lambda*ones(size(Ycell));
elseif size(lambda,1)==1
    lambda = repmat(lambda,size(Ycell,1),1);
end

if size(lambda)~=size(Ycell)
    error('lambda must have same number of dimensions as Ycell.');
end

    
end
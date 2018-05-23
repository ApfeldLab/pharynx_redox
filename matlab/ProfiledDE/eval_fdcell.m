function y = eval_fdcell(t,fd_cell,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function eval_fdcell
%
% Evaluates a cell array of functional data objects fd_cell at times t
% (possibly a cell array) with derivative k (possibly an array). Returns a 
% matrix whose columns correspond to the elements of fd_cell.  
%
% INPUTS
%
% t  -  the time points at which to evaluate the functional data objects. 
%       If given as a cell array with the same dimension as fd_cell, each
%       element of fd_cell has a unique set of evaluation points.  If t has 
%       1 column but the same number of rows as fd_cell then columns share
%       evaluation points within rows of fd_cell.  If t is a vector then it
%       is used as a common set of evaluation points for all elements of
%       fd_cell. 
%
% fd_cell - a cell array of functional data objects to be evaluated. 
%
% k  -  the order of derivative to evaluate. If k is given as  a vector,
%       but fd_cell is not a unidimensional array, k is assumed to be
%       constant within rows of fd_cell. If k is a singleton, it is common
%       to all elements of fd_cell. If it is not given, it is assumed to be
%       zero. 
%
% OUTPUT
%
% y  -  a cell array of vectors providing the evaluation of each element of
%       fd_cell at the corresponding time points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    k = 0;
end

y = cell(size(fd_cell));

t2 = t;

if length(k)==1
    k = k*ones(size(fd_cell));  
end
if ~prod(+(size(k)==size(fd_cell)))    % true if dimensions are not the same
    if min(size(k))==1 && max(size(k)~=1)
        k=reshape(k,length(k),1);
        k=repmat(k,1,size(fd_cell,2));
    else
        error('k must be the same size as the number of equations you have')
    end
end

%if t is a cell with the same number of rows as fd_cell
if iscell(t) &&  size(t,1)==size(fd_cell,1) && size(t,2)==size(fd_cell,2)
    for i = 1:numel(fd_cell)
        if ~isempty(t{i})
            y{i} = eval_fd(t{i},fd_cell{i},k(i));
        end
    end
%if t is a cell column and should be constant within rows    
elseif iscell(t) &&  size(t,1)==size(fd_cell,1) && size(t,2)==1 
    for i = 1:size(fd_cell,1)    
        for j = 1:size(fd_cell,2) 
            if ~isempty(t2{i}) 
                y{i,j} = eval_fd(t2{i},fd_cell{i,j},k(i,j));
            end
        end
    end
% if t is a vector    
elseif ~iscell(t) && min(size(t2))==1 
    for i = 1:numel(fd_cell) 
        y{i} = eval_fd(t2,fd_cell{i},k(i));
    end
else
    error('t is the wrong dimension')
end


end
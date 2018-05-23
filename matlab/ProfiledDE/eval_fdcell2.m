function y = eval_fdcell2(t,fd_cell,k)
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

y = cell2mat(eval_fdcell(t,fd_cell,k));


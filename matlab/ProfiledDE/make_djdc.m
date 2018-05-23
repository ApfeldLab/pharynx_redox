function dJ = make_djdc(fd_cell,Ycell,Tcell,wts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function djdc_pts
%
% Calculates the derivative of squared error loss with respect to the
% co-efficients in fd_cell.
%
% INPUTS:
%
% DEfd    - a cell array, each component contains a functional data 
%           object corresponding to the same component of the DE. 
%
% Ycell  -  a cell array of observations, one element corresponding to one
%           measured component. 
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% if wts is empty, the inverse standard deviation of each component
% is used as a weighting for each dimension. If lambda is a singleton, it 
% is multiplied by the average weight in each component. 
%
% OUTPUT:
%
% dJ    - the derivative of squared error loss with respect to the
%         co-efficients in fd_cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f,J] = djdc_pts(fd_cell,Ycell,Tcell,wts);

dJ = -J'*f;

end
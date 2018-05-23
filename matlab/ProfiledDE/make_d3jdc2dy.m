function d3G = make_d3jdc2dy(fd_cell,Tcell,wts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d3jdc2dy
%
% The third derivative of squared error with respect to co-efficients and
% observations when using a link function. 
%
%
% INPUTS
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% OUTPUT:
%
% d3G   -   An cell array of the same length as the co-efficients. Each 
%           element contains the third derivative matrix of the spline 
%           criterion with respect to the coefficients and observations and
%           one more co-efficient.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d3G = 0;

end
function [basis_cell Lfd_cell quadvals_cell, nbasis] = MakeBasis_cell(...
    range,norder,knots_cell,nquad,dvalue,lambda0,LFD_pen_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from function MakeBasis and spiffed up to hide a for loop which
% gives cell array results for all outputs. Makes a B-spline basis with 
% quadrature points and the Lfd, both in cell array form. Sets up the Lfd
% with a LFD_pen_order (by default equal to 1) derivative penalty suitable 
% for a first approximation to the spline coefficients.  
%
% INPUTS:
%
% range  - the range of the basis should be a cell or constant
% 
% norder - order of the basis functions can be a matrix
%
% knots_cell  - knots for the bspline basis should be a cell
%
% dvalue   -  order of derivative to store quadrature values up to.  
%
% lambda0   - smoothing parameter, if it's a single number lambda0 is kept
%             constant for all runs and components.  If lambda0 is a row 
%             vector it should be of length (number of components) and is
%             assumed constant across runs.  If lambda0 is a column vector
%             it should be of length (number of runs) and is assumed
%             constant across components with runs.  
%
% LFD_pen_order  - (optional) derivative order for the penalty to put on
%                   the basis when smoothing.  Default value is the smaller
%                   of 1 and (norder-1).
%
% OUTPUT:
%
% basis_cell    - a B-spline basis cell array with quadrature values 
%                 attached.
%
% Lfd_cell     - a Lfd object cell array.
%
% quadvals_cell -  quadrature points to put as values for the bases again a 
%                  cell array
%
% nbasis    - a cell array containing the number of nasis functions.  I
%             doubt it's really worth keeping...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<7
    LFD_pen_order=[];
end
if ~iscell(range)
    temp     = range;
    range    = cell(1,1);
    range{1} = temp;
    range    = repmat(range,size(knots_cell,1),1);% only replicate along number of runs since it should be constant along components.
end

if length(norder)==1
    norder=norder*ones(size(knots_cell));
end

%sort out lambda0
if length(lambda0)==1
    lambda0 = lambda0*ones(size(knots_cell));
elseif size(lambda0,1)==size(knots_cell,1) & size(lambda0,2)==1&length(lambda0)~=1 %i.e. same number of rows but only 1 column
    lambda0 = repmat(lambda0,1,size(knots_cell,2));
elseif size(lambda0,2)==size(knots_cell,2) & size(lambda0,1)==1 &length(lambda0)~=1 %i.e. same number of columns but only 1 row
    lambda0 = repmat(lambda0,size(knots_cell,1),1);
end
if  ~prod(+(size(lambda0)==size(knots_cell))))
    error('Wrong number of lambda0s');
end

% set up outputs
nbasis                = cell(size(knots_cell));
basisvalues           = cell(dvalue+1,1);
basis_cell            = cell(size(knots_cell));
Lfd_cell              = cell(size(knots_cell));
quadvals_cell         = cell(size(knots_cell,1),1);
make_quads_with_these = cell(size(knots_cell,1),1);  % Sometimes called BigKnots

% Gather all the knots from all the components within a run.  Use them
% later to get quad points.
for i=1:size(knots_cell,1)
    for j=1:size(knots_cell,2)
        make_quads_with_these{i} = [make_quads_with_these{i};reshape(knots_cell{i,j},length(knots_cell{i,j}),1)];
    end
end

 
for i=1:size(knots_cell,1)
    % Create simpson's rule quadrature points and values
    quadvals_cell{i} = MakeQuadPoints(make_quads_with_these{i},nquad);  
    for j=1:size(knots_cell,2)
        nbasis{i,j} = length(knots_cell{i,j}) + norder(i,j) -2;
        basis_cell{i,j} = MakeBasis(range{i},nbasis{i,j},norder(i,j),knots_cell{i,j},quadvals_cell{i},(min((dvalue+1),norder(i,j))-1));
        if isempty(LFD_pen_order)
            Lfd_cell{i,j} = fdPar(basis_cell{i,j},max(0,min(1,norder(i,j)-2)),lambda0(i));         % pts  attatched
        else
            Lfd_cell{i,j} = fdPar(basis_cell{i,j},min(LFD_pen_order,max(0,norder(i,j)-2)),lambda0(i));         % pts  attatched
        end
    end
end
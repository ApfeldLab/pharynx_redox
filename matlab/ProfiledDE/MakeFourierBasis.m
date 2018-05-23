function basis = MakeFourierBasis(range,nbasis,period,quadvals,dvalue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function MakeBasis
%
% Makes a B-spline basis with quadrature points.
%
% INPUTS:
%
% range  - the range of the basis
% 
% nbasis - number of basis points
%
% period - period of the basis function
%
% quadvals -  quadrature points to put as values for the bases
%
% dvalue   -  order of derivative to store quadrature values up to.
%
% OUTPUT:
%
% basis    - a B-spline basis with quadrature values attached.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basis = create_fourier_basis(range,nbasis,period);

basis = putquadvals(basis,quadvals);

clear values
for ivalue=1:(dvalue+1)
    basisvalues    = eval_basis(quadvals(:,1), basis, ivalue-1);
    values{ivalue} = basisvalues;
end

basis = putvalues(basis, values);

end
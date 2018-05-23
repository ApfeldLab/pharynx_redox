function [f,J] = SplineParsErr(pars,DEfd,fn,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SplineParsErr
%
% Residuals and derivatives for fitting parameters to a DE smooth. 
%
% INPUTS
%
% pars  -   initial vector of parameters
%
% DEfd  -   a cell array of functional data objects representing the smooth
%           of data to a DE.
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:   a handle to a function of the form f(t,y,p) giving the right 
%           hand side of the DE, should accept vector t, cell array 
%           function data object, fd_cell, and output cell array of
%           vectors, one for each component. 
%
%  fn.dfdp: a handle to a function giving the derivative of fn.fn with 
%           respect to the parameters. Output should be a cell array of 
%           vectors, the first dimenion giving the components of fn.fn with 
%           derivatives indexed by the second dimension. 
%
% fn_extras - object containing additional input to the right hand side 
%             of the differential equation. 
%
% OUTPUT:
%
% f  -  vector residuals of the fit, and of the penalty at quadrature
%       points.
%
% J  -  derivatives of f with respect to parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    fn_extras = [];
end


DEbasis_cell = getcellbasis(DEfd);
quadvals = getquadvals(DEbasis_cell{1});
qvals = sqrt(quadvals(:,2));

dpredns = cell2mat(eval_fdcell(quadvals(:,1),DEfd,1));

fvals = cell2mat(feval_check(fn.fn,quadvals(:,1),DEfd,pars,fn_extras));

pen = (dpredns - fvals);

pen = kron(qvals,ones(1,(size(pen,2)))).*pen;

f = reshape(pen,size(pen,1)*size(pen,2),1);

if nargout > 1

   dfvals = feval_check(fn.dfdp,quadvals(:,1),DEfd,pars,fn_extras);
   
   J = cell(size(dfvals));
   
   for i = 1:size(dfvals,1)
       for j = 1:size(dfvals,2)
        J{i,j} = qvals.*dfvals{i,j};
       end
   end
   
   J = -cell2mat(J);
   
end

end
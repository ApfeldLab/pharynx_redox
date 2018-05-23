function [smooths,forces] = linforceest(basis_cell,pbasis_cell,A,...
    whichindex,lambda,lambdap,f_lfd,Tcell,Ycell,wts,force,force_extra)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function linforceest
%
% Estimates a forcing function to make a linear differential equation
% correspond to observed data. This is done by fitting the smooth and the
% forcing functions jointly. 
%
% INPUTS:
%
% basis_cell  - a cell array of basis objectss used to represent the final
%               (forced) solutions. 
% 
% pbasis_cell - a cell array of basis objects used to represen the forcing
%               functions
%
% A           - the matrix giving the co-efficients in the linear
%               differential equations
%
% whichindex  - a vector giving the components of the system which are
%               forced by unknown functions to be represented. pbasis_cell 
%               should be the same length as whichindex. 
%
% lambda      - a vector of smoothing co-efficients; will be transformed
%               into a vector if given as a scalar. 
%
% lambdap     - a vector of smoothing co-efficients for the forcing
%               functions.
%
% f_lfd       - an Lfd object specifying the smoothing penalty to penalize
%               the forcing functions. 
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% Ycell   - observed values of the DE (matrix or cell array).
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. If left
%           empty, all weights are assumed to be equal. 
% 
% force   - (optional) a function to output a possible known forcing 
%           component, must take two arguments: a vector of times, and any
%           additional input. 
%
% force_extra - additional inputs into force.
%
% OUTPUT
%
% smooth  -  a smooth of data using the forced differential equation
%
% force   -  a functional data representation of unknown forcing functions
%            for the system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% First of all check all the intput is in good order.

    if nargin<12
        force_extra = [];
    end

    if nargin<11
        force = [];
    end
    if nargin<10
        wts = ones(size(Ycell));
    end

    if isempty(whichindex)
        whichindex = 1:length(pbasis_cell);
    end
    
    [wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);
    
    if length(lambdap)==1
        lambdap = lambdap*ones(size(whichindex));
        lambdap = reshape(lambdap,numel(lambdap),1);
    end
    
    if length(lambdap)<length(basis_cell)
        lambdap2 = zeros(length(basis_cell),1);
        lambdap2(whichindex) = lambdap;
        lambdap = lambdap2;
    end

    if iscell(Ycell)
        Ycell = cell2mat(Ycell);
        Ycell = reshape(Ycell,numel(Ycell),1);
       end
    
    if iscell(wts)
        wts = cell2mat(wts);
        wts = reshape(wts,numel(wts),1);
    end
    
    
% Now form matrices    
    
    Zmat_cell = eval_basis_cell(Tcell,basis_cell,0);

    Dphimat_cell = getvalues_cell(basis_cell,1);
    phimat_cell = getvalues_cell(basis_cell,0);

    pphimat_cell = getvalues_cell(pbasis_cell,0);
    Dpphimat_cell = getvalues_cell(pbasis_cell,f_lfd);

 % Quadrature values
    
    quadvals = getquadvals(basis_cell{1});
    n = size(quadvals,1)*size(lambda,2);
    Q = spdiags(kron(lambda',quadvals(:,2)),0,n,n);
    pQ = spdiags(kron(lambdap,quadvals(:,2)),0,n,n);
    
 % Now put it all together
        
    Zmat = mattdiag_cell(Zmat_cell,0);
    Dphimat = mattdiag_cell(Dphimat_cell,0);
    phimat = mattdiag_cell(phimat_cell,0);
    
    
    wi = zeros(length(Zmat_cell),1);
    wi(whichindex)=1;
    wi = diag(wi);
    
    pphimat = [];
    Dpphimat = [];
    k = 0;
    for i=1:length(Zmat_cell)
        if wi(i,i)
          k = k+1;
          pphimat = [pphimat kron(wi(:,i),pphimat_cell{k})];
          Dpphimat = [Dpphimat kron(wi(:,i),Dpphimat_cell{k})];
        end
    end
   
    lvec = [Zmat'*diag(wts)*Ycell; zeros(size(pphimat,2),1)];
    
    A = kron(A,speye(size(Dphimat_cell{1},1),size(Dphimat_cell{1},1)));
        
    Dmat1 = (Zmat'*diag(wts)*Zmat + ( Dphimat - A*phimat)'*Q*(Dphimat-A*phimat));
    Dmat2 = ((Dphimat-A*phimat)'*Q*pphimat);
    Dmat3 = (pphimat'*Q*(Dphimat-A*phimat));
    Dmat4 = (pphimat'*Q*pphimat+Dpphimat'*pQ*Dpphimat);
    
    Dmat = [Dmat1 Dmat2; Dmat3 Dmat4];
   
    if ~isempty(force)
       fs = force(quadvals(:,1),force_extra);    
       lvec = lvec - [(Dphimat-A*phimat)'*Q*fs; - pphimat'*Q*fs];
    end
    
% Solve and form fda objects
    
    res = Dmat\lvec;
       
    c = res(1:size(Zmat,2));
    d = res((size(Zmat,2)+1):length(res));
    
    smooths = Make_fdcell(c,basis_cell);
    forces = Make_fdcell(d,pbasis_cell);

end

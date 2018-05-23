
function [f,J] = JointSplineErr(bigcoef,DEbasis_cell,coefsizes,path,tspan,lambda,fn,dfn,dfn2)

% Re-pass parameter-coefficient matrix:

pars = bigcoef(1:coefsizes(1));
newcoefs = bigcoef((coefsizes(1)+1):length(bigcoef));

% Create the appropriate functional data object

DEfd_cell = Make_fdcell(newcoefs,DEbasis_cell);

% Matrix of errors

predns = reshape(eval_fdcell(tspan,DEfd_cell,0),size(path,1),size(path,2));
err = path - predns;
err = reshape(err,size(err,1)*size(err,2),1);

% Matrix of penalty at quadrature points

quadvals = getquadvals(DEbasis_cell{1});
qvals = sqrt(quadvals(:,2));
qpredns = reshape(eval_fdcell(quadvals(:,1),DEfd_cell,0),size(quadvals,1),size(path,2));
dpredns = reshape(eval_fdcell(quadvals(:,1),DEfd_cell,1),size(quadvals,1),size(path,2));


fvals = fn(quadvals(:,1),qpredns,pars);
dfvals = dfn(quadvals(:,1),qpredns,pars);
dfvalsp = dfn2(quadvals(:,1),qpredns,pars);


pen = (dpredns - fvals);

pen = kron(qvals,ones(1,(size(pen,2)))).*pen;

pen = reshape(pen,size(pen,1)*size(pen,2),1);

% Putting it together in one long vector

f = [err; lambda * pen];

% Now to calculate a Jacobian

if nargout > 1

    Zmat_cell = eval_basis_cell(tspan,DEbasis_cell,0);
    
    Zmat = mattdiag_cell(Zmat_cell,0);
    
% Derviatives of penalty wrt coefs:

    Dphimat_cell = getvalues_cell(DEbasis_cell,1);
    phimat_cell = getvalues_cell(DEbasis_cell,0);

    for(i = 1:size(path,2))
        Dphimat_cell{i} = Dphimat_cell{i}.*...
            kron(qvals,ones(1,size(Dphimat_cell{i},2)));
        for(j = 1:size(path,2))
            dpen2{i,j} = phimat_cell{i}.*kron(dfvals(:,i,j).*qvals,...
                ones(1,size(phimat_cell{i},2)));
        end
        
    end

    dpen1 = mattdiag_cell(Dphimat_cell,0);
    dpen2 = cell2mat(dpen2);

    
    J1 = [-Zmat; lambda*(dpen1-dpen2)];
    J2 = lambda*reshape(-dfvalsp,size(dfvalsp,1)*size(fvals,2),length(pars));
    J2 = kron( kron(ones(size(fvals,2),1),qvals), ones(1,size(J2,2)) ).*J2;
    J2 = [zeros(size(Zmat,1),length(pars)); J2];
    
    J = [J2 J1];
        
end

end

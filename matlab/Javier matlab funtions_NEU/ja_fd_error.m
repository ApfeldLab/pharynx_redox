function [MSE, MSEsum] = ja_fd_error(y, argvals, fdobj)
%PLOTFIT plots discrete data along with a functional data object for
%  fitting the data.  It is designed to be used after SMOOTH.BASIS to
%  check the fit of the data offered by the FD object.
%  Note:  As of March 2, 2009, the arguments CASENAMES, VARNAMES and
%  NFINE have been removed.
%  Arguments:
%  Y         ... the data used to generate the fit
%  ARGVALS   ... discrete argument values associated with data
%  FD        ... a functional data object for fitting the data


%  Last modified 22 June 2017
%  Javier Apfeld 

% y = dataset;
% argvals = position;
% fdobj = pos_fd;

if size(argvals,1) == 1; argvals = argvals';  end

basisobj = getbasis(fdobj);
rangeval = getbasisrange(basisobj);
nbasis   = getnbasis(basisobj);


coef  = getcoef(fdobj);
coefd = size(coef);
ndim  = length(coefd);

n    = size(y,1);
nrep = coefd(2);
casenum = 1:nrep;
if ndim < 3, nvar = 1;  else nvar = coefd(3);  end


y = reshape(y, n, nrep, nvar);

%  set up number of points at which to evaluate the curves

nfine = max([201, 10*nbasis+1]);

%  set up labels for arguments, cases and variables.

fdnames   = getnames(fdobj);
argname   = fdnames{1};
casenames = fdnames{2};
varnames  = fdnames{3};

%  compute fitted values for evalargs and fine mesh of values

yhat   = reshape(eval_fd(argvals, fdobj),[n,nrep,nvar]);
res    = y - yhat;
MSE    = squeeze(mean(res.^2))';
MSE    = reshape(MSE,nvar,nrep)';
MSEsum = sum(MSE);









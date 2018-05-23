%% Profile Estimation Experiments

%% Path Functions

fn.fn = @delayfun;             % RHS function
fn.fndfdy = @delaydfdx;        % Derivative wrt inputs
fn.fndfdp = @delaydfdp;        % Derviative wrt parameters
fn.d2fdy2 = @delayd2fdx2;      % Hessian wrt inputs
fn.d2fdydp = @delayd2fdxdp;    % Hessian wrt inputs and parameters
fn.d2fdp2 = @delayd2fdp2;
fn.d3fdydp2 = @delayd3fdxdp2;
fn.d3fdy2dp = @delayd3fdx2dp;
fn.d3fdy3 = @delayd3fdx3;


%% Observation times

times = tdat;    % Observation times
tfine = tdat;        % Times to plot solutions
range = [min(tdat),max(tdat)];

Tcell = cell(1,1);
Tcell(:) = {times};      % Which components are observed at


knots = [tdat(1)-5;  tdat(853:924); max(tdat)+5];
fbasis = create_bspline_basis([min(knots) max(knots)],length(knots)+4,6,knots);
force_par = fdPar(fbasis,1,0.01);

[force_fd,beta] = smooth_monotone(tdat,fdat,force_par);

extra .force_fd = force_fd;
extra.beta = beta;


Ycell = {xdat};

%% Other parameters

pars = [ -0.5; 0.3; 0.69];

%% Fitting parameters

lambdas = 1e5;                   % Smoothing for model-based penalty

lambda0 = 1e0;   % Smoothing for 1st-derivative penalty

nknots = [160; 0; 0];    % Number of knots to use. 

nquad = 5;       % No. between-knots quadrature points. 
norder = 3;      % Order of B-spline approximation



%% Profiling optimiation control

maxit1 = 100;      % Maximum iterations interior of profiling
maxit0 = 100;     % Maximum iterations for outer optimization
mintol = 1e-6;    % Optimization tollerance

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',maxit0,'TolFun',1e-8,'TolX',1e-10); 

% Other observed optimiation control
lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',maxit0,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 

% Optimiation control within profiling
lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',maxit1,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 



%% Setting up functional data objects

addpath ('c:\matlab\fdaM')
 
%% Setting up knots  

path = Ycell;

knots_cell = {linspace(min(times),max(times),160)};
            
%% Setting up bases

basis_cell = cell(1,length(path));      % Create cell arrays. 
Lfd_cell = cell(1,length(path));

nbasis = zeros(length(path),1);

bigknots = knots_cell{1};               % bigknots used for quadrature
nbasis(1) = length(knots_cell{1}) + norder - 2;             % points

for i = 2:length(path)
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(bigknots,nquad);   % Create simpson's rule
                                             % quadrature points and values
for i = 1:length(path)
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  % create bases
        knots_cell{i},quadvals,1);                        % with quadrature
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         % pts  attatched
end

%% Smoothing the data

disp('Smoothing the data')

DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);

figure(1)
devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(Ycell)
    plot(tfine,devals{i},'g');
    hold on;
    plot(Tcell{i},Ycell{i},'x');
end
hold off;

% Set up wts

lambda = lambdas;
wts = 1;

% Smooth

disp('Re-smoothing with model-based penalty')

newpars = pars;

coefs = getcellcoefs(DEfd);

[newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr,coefs,[],[],lsopts_other,...
    basis_cell,Ycell,Tcell,wts,lambda,fn,[],newpars,extra);

DEfd = Make_fdcell(newcoefs,basis_cell);

figure(2)             % Plot results along with exact simulation
devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(Ycell)
    plot(tfine,devals{i},'g');
    hold on;
    plot(Tcell{i},Ycell{i},'x');
end
hold off


%% Profiled Estimation

disp('Profiled Penalty Estimation')


[newpars,newDEfd_cell] = Profile_GausNewt(newpars,lsopts_out,DEfd,fn,...
    lambda,Ycell,Tcell,wts,[],lsopts_in,extra);
disp(newpars);



%% Re-smooth with Profile-Estimated Parameters

disp('A final smooth with final parameters')

figure(3)
devals = eval_fdcell(tfine,newDEfd_cell,0);
for i = 1:length(Ycell)
    plot(tfine,devals{i},'g');
    hold on;
    plot(Tcell{i},Ycell{i},'x');
end
hold off

%% Estimate Variance


d2Jdp2 = make_d2jdp2(newDEfd_cell,fn,Tcell,lambda,pars,[],wts,...
    Ycell,extra);


d2JdpdY = make_d2jdpdy(newDEfd_cell,fn,Tcell,lambda,pars,[],wts,...
    Ycell,extra);


dpdY = -d2Jdp2\d2JdpdY;

S = make_sigma(DEfd,Tcell,Ycell,0);

Cov = dpdY*S*dpdY';

%  Standard errors

StdDev = sqrt(diag(Cov));

%  Correlations

Corr = Cov./(StdDev*StdDev');

%  Display these results

disp('Approximate covariance matrix for parameters:')
disp(num2str(Cov))

disp('Approximate standard errors of parameters:')
disp(num2str(StdDev'))

disp('Approximate correlation matrix for parameters:')
disp(num2str(Corr))

%% Profile Estimation Experiments

%% Path Functions

odefn = @fhnfun1p;       % Function for ODE solver (perturbed)
odefn = @fhnfun1;         % Function for ODE solver (exact) 
fn = @fhnfun;             % RHS function
fndfdy = @fhndfdy;        % Derivative wrt inputs (Jacobian)
fndfdp = @fhndfdp;        % Derviative wrt parameters
d2fdy2 = @fhnd2fdy2;      % Hessian wrt inputs
d2fdydp = @fhnd2fdydp;    % Hessian wrt inputs and parameters
d2fdp2 = @fhnd2fdp2;      % Hessian wrt parameters.      
d3fdy3 = @fhnd3fdy3;      % Third derivative wrt inputs.
d3fdy2dp = @fhnd3fdy2dp;  % Third derivative wrt intputs, inputs and pars.
d3fdydp2 = @fhnd3fdydp2;  % Third derivative wrt inputs, pars and pars. 
                          % dimensions = time, component, input,
                          % parameters 


%odefn = @rossfun1;         % Function for ODE solver
%fn = @rossfun;             % RHS function
%fndfdy = @rossdfdy;        % Derivative wrt inputs
%fndfdp = @rossdfdp;        % Derviative wrt parameters
%d2fdy2 = @rossd2fdy2;      % Hessian wrt inputs
%d2fdydp = @rossd2fdydp;    % Hessian wrt inputs and parameters


% odefn = @ppfun1;         % Function for ODE solver
% fn = @ppfun;             % RHS function
% fndfdy = @ppdfdx;        % Derivative wrt inputs
% fndfdp = @ppdfdp;        % Derviative wrt parameters
% d2fdy2 = @ppd2fdx2;      % Hessian wrt inputs
% d2fdydp = @ppd2fdxdp;    % Hessian wrt inputs and parameters
% 
% deltatru = 0.68;
% Nitru = 80;
% bctru = 3.3;
% bbtru = 2.25;
% kctru = 4.3;
% kbtru = 15.0;
% epsilontru = 0.25;
% lambtru = 0.4;
% mtru = 0.055;

% load 'C:\Documents and Settings\gilesh\Desktop\Jiguo's Code\Profiled_PDA\Real data\realdata068.mat'
% data = data068;
% tobs = data(:,1);
% Cobs = data(:,2);
% Bobs = data(:,3);

% pk = [deltatru, Nitru, kctru, kbtru, bctru, bbtru]; 

%% Observation times

tspan = 0:0.2:20;    % Observation times
tfine = 0:0.05:20;        % Times to plot solutions
range = [0,20];

obs_pts = cell(1,2);

obs_pts{1} = 1:length(tspan);      % Which components are observed at
%obs_pts{2} = 1:length(tspan);    % which observation times. 
obs_pts{2} = [];
% obs_pts{3} = 1:length(tspan);
% obs_pts{4} = 1:length(tspan);

% obs_pts{1} = 1:401;
% obs_pts{2} = 1:401;
% obs_pts{3} = 1:401;

 y0 = [-1,1];                            % Initial conditions
% y0 = [1.13293; -1.74953; 0.02207];
% y0 = [0.5; 0.5];
% y0 = [18.6047    1.1570    0.3310    0.3310];

 extra = [];
% extra = pk;


%% Other parameters

 pars = [0.2; 0.2; 3];           % Parameters
% pars = [-0.1; -0.2; 0.3; 0];
% pars = [epsilontru; lambtru; mtru];


 sigma = 0.5;                    % Noise Level
% sigma = 0;

% jitter = 0;                     % Perturbation for starting parameters
 jitter = 0.2;

%% Fitting parameters

% lambdas = 1e2;                   % Smoothing for model-based penalty
%lambdas = 1e4;
lambdas = 1;

lambda0 = 1e0;   % Smoothing for 1st-derivative penalty

loads = [];         % Loading factor for each component (empty means use
% loads = [1 0];    % inverse of standard error). 

nknots = [401; 0; 0];    % Number of knots to use. 
%nknots = [0; 101; 0];

nquad = 5;       % No. between-knots quadrature points. 
norder = 3;      % Order of B-spline approximation

%% Replication

nreps = 50;

covest = zeros(length(pars),length(pars),nreps);

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
    'Display','iter','MaxIter',maxit1,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 



%% First create a path

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,pars);
[plot_time,plot_path] = ode45(odefn,tfine,y0,odeopts,pars);


%% Set-up observations

tspan = cell(1,size(full_path,2));
path = tspan;

for i = 1:length(obs_pts)
    tspan{i} = full_time(obs_pts{i});
    path{i} = full_path(obs_pts{i},i);
end


%% Setting up functional data objects

addpath ('c:\matlab\fdaM')
 
%% Setting up knots  

knots_cell = KnotSelection(odefn,fn,fndfdy,range,pars,y0,nknots,20,...
                [],lambda0,odeopts,lsopts_other);

            
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


%% Now do the model-based analysis.

DEfd = smoothfd_cell(path,tspan,Lfd_cell);


lambda = lambdas(1)*ones(1,length(path)); % vectorize lambda

if isempty(loads)                     % estimate loads if not given
    loads = zeros(1,length(path));
    for i = 1:length(path) 
        loads(i) = 1./sqrt(var(path{i}));
    end
    lambda = lambda.*loads;
end

coefs = getcellcoefs(DEfd);

coefs = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_other,basis_cell,path,tspan,loads,lambda,fn,...
    fndfdy,[],[],[],pars,extra);

DEfd = Make_fdcell(coefs,basis_cell);



for g = 1:nreps

    noise_path = path;                 % Add noise
    for i = 1:length(path)
        noise_path{i} = path{i} + sigma*randn(size(path{i}));
    end

    % Set up loads

    lambda = lambdas(1)*ones(1,length(noise_path)); % vectorize lambda

    if isempty(loads)                    % estimate loads if not given
        loads = zeros(1,length(noise_path));
        for i = 1:length(noise_path)
            loads(i) = 1./sqrt(var(noise_path{i}));
        end
        lambda = lambda.*loads;
    end
  
    for h = 1:length(lambdas)

        disp([g h])
        
        lambda = lambdas(h).*loads;

        [f,J] = ProfileErr(pars,DEfd,fn,fndfdy,fndfdp,d2fdy2,d2fdydp,...
            lambda,noise_path,tspan,loads,[],[],[],[],lsopts_in,extra);

    end

    S = make_sigma(DEfd,tspan,noise_path,0);

%    S = diag(kron(loads,ones(1,101)));


    Cov = inv(J'*S*J);
    
    covest(:,:,g) = Cov;
       
    save 'fhn_parest_GNsim_lambda1.mat'  covest
end



%% Example Diagnostics -- Learning the FitzHugh-Nagumo Equations


%% Path Functions

odefn = @fhnfun1;         % Function for ODE solver (exact)

fn       = @genlinfun;        % RHS function
dfdy     = @genlindfdy;       % Derivative wrt inputs (Jacobian)
dfdp     = @genlindfdp;       % Derviative wrt parameters
d2fdy2   = @genlind2fdy2;     % Hessian wrt inputs
d2fdydp  = @genlind2fdydp;    % Cross derivatives wrt inputs and parameters
d2fdp2   = @genlind2fdp2;     % Hessian wrt parameters
d3fdy2dp = @genlind3fdy2dp;   % Third derivative wrt inputs, inputs, pars 
d3fdy3   = @genlind3fdy3;     % Third derivative wrt inputs
d3fdydp2 = @genlind3fdydp2;   % Third derivative wrt inputs, pars and pars

%% Various Parameters

y0 = [-1,1];               % Initial conditions

pars0 = [0.2; 0.2; 3];     % Parameters

pars = [0; 0; 0; 0];       % Starting parameters for estimation in linear 
                           % differential equation. 

sigma = 0.25;              % Noise Level

%% Observation times

tspan = 0:0.05:20;    % Observation times

obs_pts{1} = 1:length(tspan);      % Which components are observed at
obs_pts{2} = 1:length(tspan);      % which observation times.

tfine = 0:0.05:20;    % Times to plot solutions

%% Create paths

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,pars0);
[plot_time,plot_path] = ode45(odefn,tfine,y0,odeopts,pars0);

%% Set up observations

tspan = cell(1,size(full_path,2)); 
path = tspan;

for i = 1:length(obs_pts)
    tspan{i} = full_time(obs_pts{i});
    path{i} = full_path(obs_pts{i},i);
end

% add noise

noise_path = path;                
for i = 1:length(path)
    noise_path{i} = path{i} + sigma*randn(size(path{i}));
end

% and set wts

wts = [];

if isempty(wts)                             % estimate wts if not given
    for i = 1:length(noise_path)
        wts(i) = 1./sqrt(var(noise_path{i}));
    end
end

%% Fitting parameters
 
lambda = 1000; % Smoothing for model-based penalty
lambda = lambda*wts;

lambda0 = 1;   % Smoothing for 1st-derivative penalty

nknots = 401;    % Number of knots to use.
nquad = 5;       % No. between-knots quadrature points.
norder = 6;      % Order of B-spline approximation

%% Optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-8,'TolX',1e-10); 

% Other observed optimiation control
lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','on','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 

% Optimiation control within profiling
lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 


%% Setting up functional data objects

% set up knots

range = [min(full_time),max(full_time)];  % range of observations

knots_cell = cell(size(path));            % knots for each basis
knots_cell(:) = {linspace(range(1),range(2),nknots)};

% set up bases

basis_cell = cell(1,length(path));   % Create cell arrays.
Lfd_cell = cell(1,length(path));

nbasis = zeros(length(path),1);

bigknots = knots_cell{1};            % bigknots used for quadrature points
nbasis(1) = length(knots_cell{1}) + norder - 2;             
for i = 2:length(path)
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(bigknots,nquad);   % Create simpson's rule
                                             % quadrature points and values
for i = 1:length(path)
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  % create bases
        knots_cell{i},quadvals,4);                        % with quadrature
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         % pts  attatched
end


%% Smooth the data

DEfd = smoothfd_cell(noise_path,tspan,Lfd_cell);
coefs = getcellcoefs(DEfd);

devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(plot_time,plot_path(:,i),'b','LineWidth',2);
    hold on;
    plot(tfine,devals{i},'r','LineWidth',2);
    plot(tspan{i},noise_path{i},'b.');
    hold off;
end


%% Re-smoothing with model-based penalty

% Call the Gauss-Newton solver

[wts,lambda,path] = weightslambdaspath(wts,lambda,noise_path);

[newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr,coefs,...
    [],[],lsopts_other,basis_cell,noise_path,tspan,wts,lambda,fn,...
    dfdy,[],[],[],pars);

DEfd = Make_fdcell(newcoefs,basis_cell);

% Plot results along with exact solution

devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
    plot(tspan{i},noise_path{i},'b.');
    plot(plot_time,plot_path(:,i),'c');
    hold off
end

%% Now do the profiled estimation

[newpars,DEfd] = Profile_GausNewt(pars,lsopts_out,DEfd,...
    fn,dfdy,dfdp,d2fdy2,d2fdydp,lambda,noise_path,...
    tspan,wts,[],[],[],[],lsopts_in);

disp(newpars')

devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
    plot(tspan{i},noise_path{i},'b.');
    plot(plot_time,plot_path(:,i),'c');
    hold off
end

   
%% Calculate a Sample Information Matrix

d2Jdp2 = make_d2jdp2(DEfd,fn,dfdy,dfdp,d2fdy2,d2fdydp,d2fdp2,...
    d3fdy3,d3fdy2dp,d3fdydp2,tspan,lambda,newpars,[],wts,[],[],[],noise_path)

d2JdpdY = make_d2jdpdy(DEfd,fn,dfdy,d2fdy2,d3fdy3,dfdp,d2fdydp,d3fdy2dp,...
    tspan,lambda,newpars,[],wts,[],[],[],noise_path);


dpdY = -d2Jdp2\d2JdpdY;

S = make_sigma(DEfd,tspan,noise_path,0);

Cov = dpdY*S*dpdY'


%% Analyze lack of fit

A = reshape(newpars,2,2)';

[smooths,forces] = linforceest(basis_cell,basis_cell,A,1:2,10000,...
    0.0001,2,tspan,noise_path,wts);

ss = eval_fdcell(tfine,smooths);
fs = eval_fdcell(tfine,forces);

figure(2)
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(tfine,ss{i},'r','LineWidth',2);
    hold on;
    plot(tspan{i},noise_path{i},'b.');
    plot(plot_time,plot_path(:,i),'c');
    hold off
end


figure(3)
for i = 1:2
    subplot(2,1,i)
    plot(tfine,fs{i})
end

figure(4)
k = 0;
for i = 1:2
    for j = 1:2
        k = k+1;
        subplot(2,2,k)
        plot(ss{i},fs{j},'.')
    end
end

figure(5)
plot3(ss{1},ss{2},fs{1},'.')

X = [ss{1} ss{2} ss{1}.^2 ss{2}.^2];
inv(X'*X)*X'*fs{1}
inv(X'*X)*X'*fs{2}

%% Profile Estimation Experiments - Replicated Experiments

%% RHS Functions

odefn    = @fhnfunode;    % Function for ODE solver (exact)

fn.fn       = @fhnfun;       % RHS function
fn.dfdx     = @fhndfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @fhndfdp;      % Derviative wrt parameters
fn.d2fdx2   = @fhnd2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @fhnd2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @fhnd2fdp2;    % Hessian wrt parameters.

fn.d3fdx3   = @fhnd3fdx3;    % Third derivative wrt inputs.
fn.d3fdx2dp = @fhnd3fdx2dp;  % Third derivative wrt intputs, inputs and pars.
fn.d3fdxdp2 = @fhnd3fdxdp2;  % Third derivative wrt inputs, pars and pars. 

Afn.fn = @return_sum;
Afn.dA = zeros(1,2,5);
Afn.dA(1,1,4) = 1;
Afn.dA(1,2,5) = 1;

%% Various parameters

y0 = [-1,  1;            
       1,0.75];              

   
parind = [1:5] ;                        % Which parameters hold for which
parind = repmat(parind,size(y0,1),1);% replications
parind(2,3) = 6;
disp('Parameter indices:')
disp(num2str(parind))

pars = [0.2; 0.2; 3; 0.5; 0.75; 3];           % Parameters
disp('Parameter indices:')
disp(num2str(pars'))

sigma = 0.5;                    % Noise Level

jitter = 0.2;                                    % Perturbation for initial 
startpars = pars + jitter*randn(length(pars),1); % parameter estimates
disp('Starting parameter values: ')
disp(num2str(startpars'))

%% Observation times

tspan = 0:0.05:20;    % Observation times

obs_pts = {1:401;           % Which components are observed at
           1:2:401};      % Rows represent replications. 


tfine = 0:0.05:20;    % Times to plot solutions


%% Create trajectories

odeopts = odeset('RelTol',1e-13);

for i = 1:size(y0,1)
    [full_time(:,i),full_path(:,:,i)] = ode45(odefn,tspan,y0(i,:),odeopts,pars);
    [plot_time(:,i),plot_path(:,:,i)] = ode45(odefn,tfine,y0(i,:),odeopts,pars);
end

%% Set up observations

Tcell = cell(obs_pts);
path = cell(obs_pts);

for i = 1:size(obs_pts,1)
    for j = 1:size(obs_pts,2)
        Tcell{i,j} = full_time(obs_pts{i,j},i);
        path{i,j} = pars(4)*full_path(obs_pts{i,j},1,i)+pars(5)*full_path(obs_pts{i,j},2,i);
    end
end

% add noise

Ycell = path;
for i = 1:size(path,1)
    for j = 1:size(path,2)
        Ycell{i,j} = path{i,j} + sigma*randn(size(path{i,j}));
    end
end

% and set wts

wts = [];

if isempty(wts)                         % estimate wts if not given
    wts = ones(size(path));
    for i = 1:size(Ycell,1)
        for j = 1:size(Ycell,2)
            if ~isempty(Ycell{i,j})
                wts(i) = 1./sqrt(var(path{i,j}));
            end
        end
    end
end

%% Fitting parameters

lambda  = 1000; % Smoothing for model-based penalty 
lambda  = lambda*wts;

lambda0 = 1;    % Smoothing for 1st-derivative penalty

nknots = 401;   % Number of knots to use.
nquad  = 5;     % No. between-knots quadrature points.
norder = 3;     % Order of B-spline approximation


%% Profiling optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-8,'TolX',1e-10);

% Other observed optimiation control
lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

% Optimiation control within profiling
lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


%% Setting up functional data objects

% set up knots

range = zeros(2,2);              % Range of observations
knots_cell = cell(2,2);   % Knots for each basis

for i = 1:size(knots_cell,1)
    range(i,:) = [min(full_time(:,i)),max(full_time(:,i))];  
    knots_cell(i,:) = {linspace(range(i,1),range(i,2),nknots)};
end

% set up bases

basis_cell = cell(size(knots_cell));      % Create cell arrays.
Lfd_cell = cell(size(knots_cell));

nbasis = zeros(size(knots_cell));

bigknots = cell(size(knots_cell,1),1);    % bigknots used for quadrature points
bigknots(:) = {[]};
quadvals = bigknots;

for i = 1:size(knots_cell,1)
    for j = 1:size(knots_cell,2)
        bigknots{i} = [bigknots{i} knots_cell{i,j}];
        nbasis(i,j) = length(knots_cell{i,j}) + norder -2;
    end
%    quadvals{i} = MakeQuadPoints(bigknots{i},nquad);
    knots = unique(bigknots{i});
    n = length(knots);
    
    quadvals{i} = [knots', ones(n,1)/n];    
end



for i = 1:size(knots_cell,1)        % create bases and quadrature points
    for j = 1:size(knots_cell,2)
        basis_cell{i,j} = MakeBasis(range(i,:),nbasis(i,j),norder,...  
            knots_cell{i,j},quadvals{i},1);                        
        Lfd_cell{i,j} = fdPar(basis_cell{i,j},1,lambda0);        
    end
end



%% Smooth the data

DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);
coefs = getcellcoefs(DEfd);

devals = eval_fdcell(tfine,DEfd,0);

figure(1)
for i = 1:size(DEfd,1)
    for j = 1:size(DEfd,2)
        subplot(2,2,(i-1)*2+j)
        plot(tfine,devals{i,j},'r','LineWidth',2);
        hold on
        plot(plot_time,plot_path(:,j,i),'g');
        hold off
    end
end


%% Re-smoothing with model-based penalty


% Call the Gauss-Newton solver

[newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr_rep,coefs,[],[],...
    lsopts_other,basis_cell,Ycell,Tcell,wts,lambda,fn,[],startpars,parind,Afn);

tDEfd = Make_fdcell(newcoefs,basis_cell);

% Plot results along with exact solution

devals = eval_fdcell(tfine,tDEfd,0);
figure(2)
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j)
        plot(tfine,devals{i,j},'r','LineWidth',2);
        hold on;
        plot(plot_time,plot_path(:,j,i),'c');
        hold off
    end
end



%% Perform the Profiled Estimation

[newpars,newDEfd_cell] = Profile_GausNewt_rep(startpars,lsopts_out,parind,...
    tDEfd,fn,lambda,Ycell,Tcell,wts,[],lsopts_in,Afn);

disp(['New parameter values: ',num2str(newpars')]);



%% Plot Smooth with Profile-Estimated Parameters


devals = eval_fdcell(tfine,newDEfd_cell,0);
figure(3)
for i = 1:size(path,1)
    for j = 1:size(path,2)
        subplot(size(path,1),size(path,2),(i-1)*size(path,2)+j)
        plot(tfine,devals{i,j},'r','LineWidth',2);
        hold on;
        plot(plot_time,plot_path(:,j,i),'c');
        hold off
    end
end

 


% TESTS


[f,J] = ProfileErr_rep(startpars,parind,tDEfd,fn,lambda,Ycell,Tcell,wts,[],lsopts_in,Afn);

tJ = 0*J;
eps = 1e-6;

for i = 1:length(startpars)
    disp(i)
    tpars = startpars;
    tpars(i) = startpars(i) + eps;
    tf = ProfileErr_rep(tpars,parind,tDEfd,fn,lambda,Ycell,Tcell,wts,[],lsopts_in,Afn);
    tJ(:,i) = (tf-f)/eps;
end
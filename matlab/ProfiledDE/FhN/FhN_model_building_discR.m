%% Example Diagnostics -- Learning the FitzHugh-Nagumo Equations
%
% This page provides a demonstration of the use of forcing function
% diagnostic tools for model building in systems of differential equations.
% 
% We will use the FitzHugh-Nagumo Equations:
% 
% $$\dot{x} = c \left(x - \frac{x^3}{3} + y \right)$$
%
% $$\dot{y} = -\frac{1}{c} \left( x - a + by \right)$$
% 
% as an example, but we will will act on the basis of not knowing these
% equations. Therefore, we will start out with an autonomous linear system:
%
% $$\dot{x} = a_{11}x + a_{12}y$$
%
% $$\dot{y} = a_{21}x + a_{22}y$$
%
% and then use diagnostic tools to uncover the missing cubic term. 
%
% The format of this demonstration follows that detailed in FhNEx.html and
% commentary will therefore be restricted to diagnostic terms. 

%% RHS Functions
%
% Since we are using a linear function to begin with, we make use of the
% |genlin| set of functions (although |fhnfunode| will be used to produce
% data). See genlinEx.html for a full description. 

odefn       = @fhnfunode;        % Function for ODE solver (exact)

fn.fn       = @polyfun;        % RHS function
fn.dfdx     = @polydfdx;       % Derivative wrt inputs (Jacobian)
fn.dfdp     = @polydfdp;       % Derviative wrt parameters
fn.d2fdx2   = @polyd2fdx2;     % Hessian wrt inputs
fn.d2fdxdp  = @polyd2fdxdp;    % Cross derivatives wrt inputs and parameters
fn.d2fdp2   = @polyd2fdp2;     % Hessian wrt parameters
fn.d3fdx2dp = @polyd3fdx2dp;   % Third derivative wrt inputs, inputs, pars 
fn.d3fdx3   = @polyd3fdx3;     % Third derivative wrt inputs
fn.d3fdxdp2 = @polyd3fdxdp2;   % Third derivative wrt inputs, pars and pars


%% Various Parameters

y0 = [-1,1];               % Initial conditions

pars0 = [0.2; 0.2; 3];     % Parameters for the FitzHugh-Nagumo equations

sigma = 0.25;              % Noise Level

%% Observation times

tspan = 0:0.05:20;    % Observation times

obs_pts{1} = 1:length(tspan);      % Which components are observed at

tfine = 0:0.05:20;    % Times to plot solutions

%% Create trajectories

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,pars0);
[plot_time,plot_path] = ode45(odefn,tfine,y0,odeopts,pars0);


%% Set up observations

Tcell = cell(size(obs_pts)); 
path = Tcell;

for i = 1:length(obs_pts)
    Tcell{i} = full_time(obs_pts{i});
    path{i} = full_path(obs_pts{i},i);
end

% add noise

Ycell = path;                
for i = 1:length(obs_pts)
    Ycell{i} = path{i} + sigma*randn(size(path{i}));
end

% and set wts

wts = [];

if isempty(wts)                             % estimate wts if not given
    for i = 1:length(Ycell)
        if  ~isempty(Ycell{i}) 
            wts(i) = 1./sqrt(var(Ycell{i}));
        else
            wts(i) = 1;
        end
    end
end

%% Fitting parameters
 
lambda = 1000; % Smoothing for model-based penalty
lambda = lambda*wts;

lambda0 = 0.1;   % Smoothing for 1st-derivative penalty

nknots = 401;    % Number of knots to use.
nquad = 5;       % No. between-knots quadrature points.
norder = 6;      % Order of B-spline approximation

%% Optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',20,'TolFun',1e-8,'TolX',1e-10); 

% Other observed optimiation control
lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','on','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 

% Optimiation control within profiling
lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',20,'TolFun',1e-14,'TolX',1e-14,...
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
    Lfd_cell{i} = fdPar(basis_cell{i},3,lambda0);         % pts  attatched
end


%% Smooth the data
%
% Since we are smoothing with a Linear Differential Equation
%
% $$\dot{\mathbf{x}} = A \mathbf{x}$$
% 
% with initial parameters |A=0|, the model based penaly is equivalent to a
% first derivative penalty, and we will only use this. 

DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);
coefs = getcellcoefs(DEfd);

figure(1)
axes('fontsize',14)
devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(plot_time,plot_path(:,i),'b','LineWidth',2);
    hold on;
    plot(tfine,devals{i},'r','LineWidth',2);
    plot(Tcell{i},Ycell{i},'b.');
    hold off;
    if i==1
        ylabel('\fontsize{20} V','rotation',0)
%        title(['\fontsize{13} Raw data (.), ', ...
%               'smoothing solution (r-) and true path (g-)'])
        axis([0 20 -4 4])
    else
        axis([0 20 -2 2])
        xlabel('\fontsize{20} t')
        ylabel('\fontsize{20} R','rotation',0)
    end
end


ddevals = eval_fdcell(tfine,DEfd,1);

X = [devals{1} devals{1}.^2 devals{1}.^3 ones(size(Ycell{1}))];

startpars = (X'*X)\X'*ddevals{1};



%% Now do the profiled estimation

fn.fn       = @forcingfun;        % RHS function

fn.dfdx     = @forcingdfdx;       % Derivative wrt inputs (Jacobian)
fn.dfdp     = @forcingdfdp;       % Derviative wrt parameters

fn.d2fdx2   = @forcingd2fdx2;     % Hessian wrt inputs
fn.d2fdxdp  = @forcingd2fdxdp;    % Cross derivatives wrt inputs and parameters
fn.d2fdp2   = @forcingd2fdp2;     % Hessian wrt parameters

fn.d3fdx2dp = @forcingd3fdx2dp;   % Third derivative wrt inputs, inputs, pars 
fn.d3fdx3   = @forcingd3fdx3;     % Third derivative wrt inputs
fn.d3fdxdp2 = @forcingd3fdxdp2;   % Third derivative wrt inputs, pars and pars

%% Analyze lack of fit - estimate forcing functions. 

fn_extras.fn      = @polyfun;    % Original function
fn_extras.dfdx    = @polydfdx;   % First derivative
fn_extras.d2fdx2  = @polyd2fdx2; % Second derivative
fn_extras.d3fdx3  = @polyd3fdx3; % Third derivative

fn_extras.extras  = [];        % Original information to fn_extras.fn. 

% We also need to add parameters:

fn_extras.pars = startpars;

% We also need to specify a basis for the forcing components and a vector
% indicating which components are to be forced. 

fn_extras.basisp = basis_cell(1);   % Cell array of basis functions.
fn_extras.which  = 1;             % Force the first component of the 
                                  % system only. 

pen   = @forcingpen;   % Penalty
dpen  = @forcingdpen;  % Derivative with respect to forcing co-efficients
d2pen = @forcingd2pen; % Second derivative

% These penalty functions also require extra arguments in the form of a
% struct to specify the basis, the degree of smoothing and the smoothing
% penalty:

pen_extras.basis = fn_extras.basisp; % Same basis as fn_extras. 
pen_extras.deg    = 2;           % Penalize the second derivative
pen_extras.lambda = 0.001;      % Smoothing parameter.      


% New starting parameters

nstartpars = zeros(getnbasis(fn_extras.basisp{1}),1);


[fcoefs,nDEfd] = Profile_GausNewt(nstartpars,lsopts_out,DEfd,fn,lambda,...
    Ycell,Tcell,wts,[],lsopts_in,fn_extras,pen,dpen,pen_extras);

% DEfd now takes the form of a smooth to the data and we can see that it
% fits better:

figure(2)
ss = eval_fdcell(tfine,nDEfd,0);
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(tfine,ss{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    plot(plot_time,plot_path(:,i),'c');
    hold off
end


%% Examine Forcing Functions
%
% First we constitute the approximated forcing function:

fd_approx = fd(fcoefs,basis_cell{1});

fs = eval_fd(tfine,fd_approx);

figure(3)
plot(tfine,fs);


figure(4)
axes('fontsize',14);
plot(ss{1},fs,'linewidth',2);
xlabel('\fontsize{20} V')
ylabel('\fontsize{20} g_V','rotation',0)

%% Calculate a Sample Information Matrix

d2Jdp2 = make_d2jdp2(DEfd,fn,Ycell,Tcell,lambda,fcoefs,[],wts,...
    fn_extras,d2pen,pen_extras);

d2JdpdY = make_d2jdpdy(DEfd,fn,Ycell,Tcell,lambda,fcoefs,[],wts,fn_extras);

dpdY = -d2Jdp2\d2JdpdY;

S = make_sigma(DEfd,Tcell,Ycell,0);

Cov = dpdY*S*dpdY';



%% Calculate Hotelling Distance from Zero
%
% Look at the distance of |fcoefd| from 0 with respect to the metric
% defined by |Cov|. This is a heuristic Wald test for the goodness of fit of the 
% original equations. 

disp(['Goodness of fit = ',num2str(fcoefs'*inv(Cov)*fcoefs)])


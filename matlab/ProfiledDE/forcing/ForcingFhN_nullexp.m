%% Estimating Forcing Functions for the FitzHugh-Nagumo Equations
%
% This page provides a demonstration of the use of forcing functiona
% perturbed set of FitzHugh-Nagumo Equations:
% 
% $$\dot{x} = c \left(x - \frac{x^3}{3} + y \right) + g(t)$$
%
% $$\dot{y} = -\frac{1}{c} \left( x - a + by \right)$$
% 
% Here we will assume knowledge of $a$, $b$ and $c$, and estimate $g(t)$
% from data. 
%
% The format of this demonstration follows that detailed in FhNEx.html and
% commentary will therefore be restricted to terms unique to forcing 
% function estimation.  

%% RHS Functions
%
% Since we are using a linear function to begin with, we make use of the
% |forcing| set of functions (although |fhnfunodep| will be used to produce
% data). 

odefn       = @fhnfunode;        % Function for ODE solver (exact)

fn.fn       = @forcingfun;        % RHS function

fn.dfdx     = @forcingdfdx;       % Derivative wrt inputs (Jacobian)
fn.dfdp     = @forcingdfdp;       % Derviative wrt parameters

fn.d2fdx2   = @forcingd2fdx2;     % Hessian wrt inputs
fn.d2fdxdp  = @forcingd2fdxdp;    % Cross derivatives wrt inputs and parameters
fn.d2fdp2   = @forcingd2fdp2;     % Hessian wrt parameters

fn.d3fdx2dp = @forcingd3fdx2dp;   % Third derivative wrt inputs, inputs, pars 
fn.d3fdx3   = @forcingd3fdx3;     % Third derivative wrt inputs
fn.d3fdxdp2 = @forcingd3fdxdp2;   % Third derivative wrt inputs, pars and pars

%% Observation times

tspan = 0:0.05:20;    % Observation times

obs_pts{1} = 1:length(tspan);      % Which components are observed at
obs_pts{2} = 1:length(tspan);      % which observation times.

tfine = 0:0.05:20;    % Times to plot solutions

%% Various Parameters

y0 = [-1,1];               % Initial conditions
pars0 = [0.2; 0.2; 3];     % Parameters for the FitzHugh-Nagumo equations

% Set up a perturbation functional data object

basis_obj = create_bspline_basis([0 20],42,3,0:0.5:20);
quadvals = MakeQuadPoints(0:0.5:20,21);      % We will need to use
basis_obj = putquadvals(basis_obj,quadvals); % quadrature points later


% Finally decide on a noise level:

sigma = 0.25;              % Noise Level

%% Extra Information for the System:
%
% In particular, we need to specify the original ODEs and their derifatives
% with respect to components (in this case the FitzHugh-Nagumo Equations)
% and also the option of adding further information for the original
% system. 

fn_extras.fn      = @fhnfun;    % Original function
fn_extras.dfdx    = @fhndfdx;   % First derivative
fn_extras.d2fdx2  = @fhnd2fdx2; % Second derivative
fn_extras.d3fdx3  = @fhnd3fdx3; % Third derivative

fn_extras.extras  = [];        % Original information to fn_extras.fn. 

% We also need to add parameters:

fn_extras.pars = pars0;

% We also need to specify a basis for the forcing components and a vector
% indicating which components are to be forced. 

fn_extras.basisp = {basis_obj};   % Cell array of basis functions.
fn_extras.which  = 1;             % Force the first component of the 
                                  % system only. 

% Note that although we are only forcing one component here, the forcing 
% basis is still represented as a cell array and we could equally have 
% estimated a number of forcing functions.         
                                  
%% Penalties on Forcing Functions
%
% We can also place roughness penalties on forcing functions, these will
% then occur as inputs into |Profile_GausNewt|. 

pen   = @forcingpen;   % Penalty
dpen  = @forcingdpen;  % Derivative with respect to forcing co-efficients
d2pen = @forcingd2pen; % Second derivative

% These penalty functions also require extra arguments in the form of a
% struct to specify the basis, the degree of smoothing and the smoothing
% penalty:

pen_extras. basis = {basis_obj}; % Same basis as fn_extras. 
pen_extras.deg    = 2;           % Penalize the second derivative
pen_extras.lambda = 0.0001;      % Smoothing parameter.                            

%% Initial Forcing Estimates
%
% We start off assuming the forcing function is zero. Since the
% coefficients of a basis expansion for it occur as parameters in the
% profiled estimation scheme we set them as:

startpars = zeros(getnbasis(basis_obj),1);

%% Create trajectories

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,pars0);
[plot_time,plot_path] = ode45(odefn,tfine,y0,odeopts,pars0);

%% Set up observations

Tcell = cell(1,size(full_path,2)); 
path = Tcell;

for i = 1:length(obs_pts)
    Tcell{i} = full_time(obs_pts{i});
    path{i} = full_path(obs_pts{i},i);
end



% and set wts

wts = [];

if isempty(wts)                             % estimate wts if not given
    for i = 1:length(path)
        wts(i) = 1./sqrt(var(path{i}));
    end
end

%% Fitting parameters
 
lambdas = exp(-2:4);          % Smoothing for model-based penalty


lambda0 = 1;   % Smoothing for 1st-derivative penalty

nknots = 401;    % Number of knots to use.
nquad = 5;       % No. between-knots quadrature points.
norder = 6;      % Order of B-spline approximation


nrep = 50;
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

DEfd = smoothfd_cell(path,Tcell,Lfd_cell);
coefs = getcellcoefs(DEfd);


%% Re-Smooth with a Model-Based Penalty

[newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_other,basis_cell,path,Tcell,wts,lambdas(1),fn,[],startpars,...
    fn_extras);

tDEfd = Make_fdcell(newcoefs,basis_cell);

htds = zeros(nrep,length(lambdas));

for g = 1:nrep
    
    Ycell = path;
    for i = 1:length(path)
        Ycell{i} = path{i} + sigma*randn(size(path{i}));
    end
    
    for h = 1:length(lambdas)

        lambda = lambdas(h);


        [fcoefs,DEfd] = Profile_GausNewt(startpars,lsopts_out,tDEfd,fn,lambda,...
            Ycell,Tcell,wts,[],lsopts_in,fn_extras,pen,dpen,pen_extras);


        %% Calculate a Sample Information Matrix

        d2Jdp2 = make_d2jdp2(DEfd,fn,Ycell,Tcell,lambda,fcoefs,[],wts,...
            fn_extras,d2pen,pen_extras);

        d2JdpdY = make_d2jdpdy(DEfd,fn,Ycell,Tcell,lambda,fcoefs,[],wts,fn_extras);

        dpdY = -d2Jdp2\d2JdpdY;

        S = make_sigma(DEfd,Tcell,Ycell,0);

        Cov = dpdY*S*dpdY';

        %% Calculate Hotelling Distance from Zero

        ht2 = fcoefs'*inv(Cov)*fcoefs;
        
        disp(['Goodness of fit measures: ',num2str([g h ht2])])
        
        htds(g,h) = ht2;
    end
end

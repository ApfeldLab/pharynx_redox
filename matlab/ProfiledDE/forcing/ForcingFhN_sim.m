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

odefn       = @fhnfunodep;        % Function for ODE solver (exact)

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

pbasis = create_bspline_basis([0 20],3,0,[0 7 14 20]);

% Finally decide on a noise level:

sigma = 0.25;              % Noise Level

%% Extra Information for the System:

fn_extras.fn      = @fhnfun;    % Original function
fn_extras.dfdx    = @fhndfdx;   % First derivative
fn_extras.d2fdx2  = @fhnd2fdx2; % Second derivative
fn_extras.d3fdx3  = @fhnd3fdx3; % Third derivative

fn_extras.extras  = [];        % Original information to fn_extras.fn. 

fn_extras.pars = pars0;

% We also need to specify a basis for the forcing components and a vector
% indicating which components are to be forced. 

fn_extras.basisp = {basis_obj};   % Cell array of basis functions.
fn_extras.which  = 1;             % Force the first component of the 
                                  % system only. 

%% Penalties on Forcing Functions

pen   = @forcingpen;   % Penalty
dpen  = @forcingdpen;  % Derivative with respect to forcing co-efficients
d2pen = @forcingd2pen; % Second derivative

pen_extras. basis = {basis_obj}; % Same basis as fn_extras. 
pen_extras.deg    = 2;           % Penalize the second derivative

%% Initial Forcing Estimates


startpars = zeros(getnbasis(basis_obj),1);



%% Fitting parameters
 
lambdas = 10.^8;
lambdap = 10.^(-4:2:0);
levels = 0;%:0.5:2;
nrep = 50;

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


lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun); 


%% Setting up functional data objects

% set up knots

range = [0,20];  % range of observations

knots_cell = cell(1,2);            % knots for each basis
knots_cell(:) = {linspace(range(1),range(2),nknots)};

% set up bases

basis_cell = cell(1,2);   % Create cell arrays.
Lfd_cell = cell(1,2);

nbasis = zeros(2,1);

bigknots = knots_cell{1};            % bigknots used for quadrature points
nbasis(1) = length(knots_cell{1}) + norder - 2;             
for i = 2:2
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(bigknots,nquad);   % Create simpson's rule
                                             % quadrature points and values
for i = 1:2
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  % create bases
        knots_cell{i},quadvals,4);                        % with quadrature
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         % pts  attatched
end


wald = zeros(length(levels),length(lambdas),length(lambdap),nrep);

sse1 = zeros(length(levels),length(lambdas),nrep);
sse2 = zeros(length(levels),length(lambdas),length(lambdap),nrep);
lr = zeros(length(levels),length(lambdas),length(lambdap),nrep);
dhYdY = zeros(length(startpars),length(startpars),length(levels),length(lambdas),length(lambdap),nrep);
dpdY = zeros(length(startpars),length(startpars),length(levels),length(lambdas),length(lambdap),nrep);


for h = 1:length(levels)

    %% Create trajectories

    fd_obj = fd([0 levels(h) 0]',pbasis);

    odeopts = odeset('RelTol',1e-13);
    [full_time,full_path] = ode45(odefn,tspan,y0,odeopts,pars0,fd_obj);
    [plot_time,plot_path] = ode45(odefn,tfine,y0,odeopts,pars0,fd_obj);

    %% Set up observations

    Tcell = cell(1,size(full_path,2));
    path = Tcell;

    for g = 1:nrep
        
        for i = 1:length(obs_pts)
            Tcell{i} = full_time(obs_pts{i});
            path{i} = full_path(obs_pts{i},i);
        end

        % add noise

        Ycell = path;
        for i = 1:length(path)
            Ycell{i} = path{i} + sigma*randn(size(path{i}));
        end

        % and set wts

        wts = [1 1];

        DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);
        coefs = getcellcoefs(DEfd);
        
        %% Now do the profiled estimation
        
        for i = 1:length(lambdas)

            [newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr,coefs,[],[],...
                lsopts_other,basis_cell,Ycell,Tcell,wts,lambdas(i),fn,[],...
                startpars,fn_extras);

            tDEfd = Make_fdcell(newcoefs,basis_cell);
            
            sse1(h,i,g) = sum( sum( (cell2mat(Ycell) - cell2mat(eval_fdcell(Tcell,tDEfd))).^2 ));
            
            for j = 1:length(lambdap)

                disp([h i j g])
                
                pen_extras.lambda = lambdap(j);
                
                [fcoefs,DEfd] = Profile_GausNewt(startpars,lsopts_out,...
                    DEfd,fn,lambdas(i),Ycell,Tcell,wts,[],lsopts_in,...
                    fn_extras,pen,dpen,pen_extras);

                %% Calculate a Sample Information Matrix

                [dhYdY,dpdY,df(:,h,i,j,g)] = est_df(DEfd,fn,Ycell,Tcell,...
                    wts,lambdas(i),fcoefs,[],fn_extras);
                
                Cov = dpdY*dpdY';

                wald(h,i,j,g) = fcoefs'*inv(Cov)*fcoefs;
                
                sse2(h,i,j,g) = sum( sum( (cell2mat(Ycell) - cell2mat(eval_fdcell(Tcell,DEfd))).^2 ));
                
                lr(h,i,j,g) = sse1(h,i,g) - sse2(h,i,j,g);
                
            end
        end
    end
    
    save '../../experiments/forcingfhn_051206.mat'
end


mwald = mean(wald,4);
mlr = mean(lr,4);
msse2 = mean(sse2,4);
msse1 = mean(sse1,3);


wald05 = quantile(squeeze(wald(1,:,:,:)),0.95,2);
lr05 = quantile(squeeze(lr(1,:,:,:)),0.95,2);

pwald = zeros(length(levels),length(lambdas),length(lambdap));
plr = zeros(length(levels),length(lambdas),length(lambdap));

for i = 1:length(lambdas)
    for j = 1:length(lambdap)
        pwald(:,i,j) = sum(squeeze(wald(:,i,j,:))>wald05(i,j),2);
        plr(:,i,j) = sum(squeeze(lr(:,i,j,:))>lr05(i,j),2);
    end
end
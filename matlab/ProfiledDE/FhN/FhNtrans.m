%% Demonstration of Profiled Estimation of Differential Equations

odefn       = @fhnfunode;      % Function for ODE solver (exact)
fn.fn       = @fhnfun;       % RHS function

% Now derivatives of the function with respect to system components and
% parameters 

fn.dfdx     = @fhndfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @fhndfdp;      % Derviative wrt parameters

% Now we need functions to compute all three sets of second derivatives:

fn.d2fdx2   = @fhnd2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @fhnd2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @fhnd2fdp2;    % Hessian wrt parameters.    

% Finally, if we want to have variance estimates, various third derivatives
% are needed:

fn.d3fdx3   = @fhnd3fdx3;    % Third derivative wrt inputs.
fn.d3fdx2dp = @fhnd3fdx2dp;  % Third derivative wrt intputs, inputs and pars.
fn.d3fdxdp2 = @fhnd3fdxdp2;  % Third derivative wrt inputs, pars and pars. 
                             % dimensions = time, component, input,
                             % parameters
                         

% Now we'll define a linear transformation of the data

Afn.fn = @return_sum;
Afn.dA = zeros(1,2,5);
Afn.dA(1,1,4) = 1;
Afn.dA(1,2,5) = 1;
                             
%% Various Parameters

% In order to specify a solution to a set of differential equations, we
% need to know initial conditions:

y0 = [-1,1];                      

% and parameters

pars = [0.2; 0.2; 3; 0.5; 0.75];          
disp(['Parameter values: ',num2str(pars')])

% We also need to specify the variance of observational noise:

sigma = 0.1;
jitter = 0.2;                   

startpars = pars + jitter*randn(length(pars),1);
disp(['Initial parameter values: ',num2str(startpars')])

%% Observation times

tspan = 0:0.05:20;

obs_pts{1} = 1:length(tspan);       

tfine = 0:0.05:20;     

%% Calculate trajectories

odeopts = odeset('RelTol',1e-13);

[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,pars);
[plot_time,plot_path] = ode45(odefn,tfine,y0,odeopts,pars);

%% Set up observations

% We start by defining cell arrays:

Tcell = cell(1,length(obs_pts));
path_cell = Tcell;

% We take the data from the solution of the differential equation and put
% it into the appropriate component. 


Tcell{1} = full_time(obs_pts{1});
path_cell{1} = pars(4)*full_path(:,1) + pars(5)*full_path(:,2);

% Finally, we add random observational noise to the 'path' variable. 

Ycell = path_cell;                
for i = 1:length(path_cell)
    Ycell{i} = path_cell{i} + sigma*randn(size(path_cell{i}));
end

% To observe the raw data and the true path, we can plot

figure(1)
for i = 1:length(Ycell)
    subplot(length(Ycell),1,i)
    plot(tspan,path_cell{i},'r')
    hold on
    plot(Tcell{i},Ycell{i},'b.')
    hold off
    if i==1
        ylabel('\fontsize{13} V')
        title('\fontsize{13} Raw data (.) and true path (-)')
    else
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} R')
    end
end

wts = [];  
                 
for i = 1:length(Ycell)
    if ~isempty(Ycell{i})
        wts(i) = 1./sqrt(var(path_cell{i}));
    else
        wts(i) = 1;
    end
end

%% Fitting parameters

lambda = 1e4;      
lambda = lambda *ones(1,2);


lambda0 = 1;        

nknots = 401;    
norder = 3;
nquad = 5;     

%% Profiling optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-6,'TolX',1e-10);

lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


%% Setting up Functional Data Objects

range = [min(full_time),max(full_time)];

knots_cell = cell(1,2);
knots_cell(:) = {linspace(range(1),range(2),401)};

basis_cell = cell(1,2); 

Lfd_cell = cell(1,2);

nbasis = zeros(2,1);

bigknots = knots_cell{1};               
nbasis(1) = length(knots_cell{1}) + norder - 2;          

for i = 2:length(knots_cell)
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

%quadvals = MakeQuadPoints(bigknots,nquad);   

knots = unique(bigknots);
n = length(knots);

quadvals = [knots', ones(n,1)/n];

for i = 1:length(knots_cell)
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,knots_cell{i},quadvals,1);                        
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         
end


%% Smoothing with model-based penalty

DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);

figure(3)
devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(devals)
    subplot(length(devals),1,i);
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on
    plot(plot_time,plot_path(:,i),'g');
    hold off
end


coefs = getcellcoefs(DEfd);

[newcoefs,resnorm2] = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_out,basis_cell,Ycell,Tcell,wts,lambda,fn,[],startpars,Afn);

DEfd = Make_fdcell(newcoefs,basis_cell);

% plot results along with exact solution

figure(3)
devals = eval_fdcell(tfine,DEfd,0);
for i = 1:length(devals)
    subplot(length(devals),1,i);
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on
    plot(plot_time,plot_path(:,i),'g');
    hold off
end

figure(4)
plot(tfine,startpars(4)*devals{1}+startpars(5)*devals{2},'r')
hold on
plot(Tcell{1},Ycell{1},'.')
plot(tspan,path_cell{1},'g')
hold off




%% Perform the Profiled Estimation

[newpars,newDEfd_cell] = Profile_GausNewt(startpars,lsopts_out,DEfd,fn,...
    lambda,Ycell,Tcell,wts,[],lsopts_in,Afn,[],[]);

disp(['New parameter estimates: ',num2str(newpars')]);


% plot smooth with profile-estimated parameters


figure(3)
devals = eval_fdcell(tfine,newDEfd_cell,0);
for i = 1:length(devals)
    subplot(length(devals),1,i);
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on
    plot(plot_time,plot_path(:,i),'g');
    hold off
end

figure(4)
plot(tfine,newpars(4)*devals{1}+newpars(5)*devals{2},'r')
hold on
plot(Tcell{1},Ycell{1},'.')
plot(tspan,path_cell{1},'g')
hold off


% Squared error for truth

F = ProfileErr(pars,DEfd,fn,lambda,Ycell,Tcell,wts,[],lsopts_other,Afn);

disp(F'*F);


F2 = ProfileErr(newpars,newDEfd_cell,fn,lambda,Ycell,Tcell,wts,[],lsopts_other,Afn);

disp(F2'*F2);

% Checking stuff

[f,J] = ProfileErr(startpars,DEfd,fn,lambda,Ycell,Tcell,wts,[],lsopts_in,Afn);

tJ = 0*J;
eps = 1e-6;

for i = 1:length(startpars)
    tpars = startpars;
    tpars(i) = startpars(i) + eps;
    tf = ProfileErr(tpars,DEfd,fn,lambda,Ycell,Tcell,wts,[],lsopts_in,Afn);
    tJ(:,i) = (tf-f)/eps;
end


d2Gdc2 = make_d2gdc2(DEfd,fn,Tcell,wts,lambda,startpars,alg,Afn);

d2Gdcdp = make_d2gdcdp(DEfd,fn,lambda,startpars,alg,Afn,Ycell,Tcell,wts);
    
dcdp = -d2Gdc2\d2Gdcdp;

coefs = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_out,basis_cell,Ycell,Tcell,wts,lambda,fn,[],startpars,Afn);


d2G = 0*d2Gdcdp;


A = Afn.fn(startpars);
[f,Zmat] = djdc_pts(DEfd,Ycell,Tcell,wts,A);
dgdc = -Zmat'*f;

for i = 1:length(startpars)
    tpars = startpars;
    tpars(i) = startpars(i) + eps;    
    tA = Afn.fn(tpars);
    [terr,tzmat] = djdc_pts(DEfd,Ycell,Tcell,wts,tA);
    
    tdgdc = -tzmat'*terr;
    d2G(:,i) = (tdgdc-dgdc)/eps;
end

d2H = make_d2SSEdcdp(DEfd,pars,Afn,Ycell,Tcell,wts);


tG = zeros(length(coefs),length(pars));

for i = 1:length(startpars)
    tpars = startpars;
    tpars(i) = startpars(i) + eps;
    tcoefs = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_out,basis_cell,Ycell,Tcell,wts,lambda,fn,[],tpars,Afn);

    tG(:,i) = (tcoefs-coefs)/eps;
end




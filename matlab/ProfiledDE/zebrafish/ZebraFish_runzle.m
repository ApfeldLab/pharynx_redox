% Model-Building for spiking rate of a Zebra Fish/

% RHS Functions

fn.fn       = @runzlefun;       % RHS function

% Now derivatives of the function with respect to system components and
% parameters 

fn.dfdx     = @runzledfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @runzledfdp;      % Derviative wrt parameters

% Now we need functions to compute all three sets of second derivatives:

fn.d2fdx2   = @runzled2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @runzled2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @runzled2fdp2;    % Hessian wrt parameters.    


startpars = [4.9597
             0.2629
            -0.0092
            -0.0006
           -11.6474
            -0.4643
             0.0003
            -0.0075];
    

%% Load the data

load 'test03.asc'

Ycell = cell(1,2);
Tcell = cell(1,2);

Ycell{1} = test03(2000:3000,1);
Tcell{1} = 0:1000;

tfine = 0:1000;

% To observe the raw data and the true path, we can plot

figure(1)
axes('fontsize',14)
plot(Tcell{1},Ycell{1},'.')
xlabel('\fontsize{20} t')
ylabel('\fontsize{20} V','rotation',0)

wts = [1 1];


%% Fitting parameters

lambda = 1e6;      
lambda = lambda * wts;

lambda0 = 1;        

nknots = 501;    
norder = 6;
nquad = 5;     

%% Profiling optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',100,'TolFun',1e-8,'TolX',1e-10);

lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',200,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',100,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


%% Setting up Functional Data Objects

range = [0 1000];

nbasis = nknots + norder - 2;        
knots = linspace(range(1),range(2),nknots);

quadvals = MakeQuadPoints(knots,nquad);   

basis_obj = MakeBasis(range,nbasis,norder,knots,quadvals,3);
Lfd_obj = fdPar(basis_obj,3,lambda0);

Lfd_cell = {Lfd_obj Lfd_obj};


%% Smooth the data

DEfd1 = smoothfd_cell(Ycell,Tcell,Lfd_cell);

% plot the smooth

figure(2)
axes('fontsize',14)
devals = eval_fdcell(tfine,DEfd1,0);
plot(tfine,devals{1},'r','LineWidth',2);
hold on;
plot(Tcell{1},Ycell{1},'b.');
hold off;

xlabel('\fontsize{13} t')
ylabel('\fontsize{13} V','rotation',0)


% Use co-efficients as a starting point

coefs = getcellcoefs(DEfd1);

% We also want to estimate parameters from this smooth

ddevals = eval_fdcell(tfine,DEfd1,1);

figure(3)
axes('fontsize',14)
plot(devals{1},ddevals{1},'linewidth',2);
xlabel('\fontsize{20} V')
ylabel('\fontsize{20} dV/dt','rotation',90)

d2devals = eval_fdcell(tfine,DEfd1,2);
plot3(devals{1},ddevals{1},d2devals{1});

%% Model-based smooth:

basis_cell = getcellbasis(DEfd1);

startcoefs2 = lsqnonlin(@SplineCoefErr_DEfit,getcoef(DEfd1{2}),[],[],...
    lsopts_out,DEfd1,2,fn,startpars,[]);

DEfd2 = update_fdcell(startcoefs2,2,DEfd1);


newcoefs = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_out,basis_cell,Ycell,Tcell,wts,lambda,fn,[],startpars);

tDEfd = Make_fdcell(newcoefs,basis_cell);

% plot results along with exact solution

figure(3)
devals = eval_fdcell(tfine,tDEfd,0);
axes('fontsize',14)
for i = 1:length(Ycell)
    subplot(length(Ycell),1,i);
    plot(Tcell{i},Ycell{i},'b.');
    hold on
    plot(tfine,devals{i},'r');
    hold off
end


ddevals = eval_fdcell(tfine,tDEfd,1);
fdevals = runzlefun(tfine,tDEfd,startpars);

figure(4)
subplot(2,1,1)
plot(ddevals{1}-fdevals{1})
subplot(2,1,2)
plot(ddevals{2}-fdevals{2})


figure(5)
for i = 1:2 
   for j = 1:2
      subplot(2,2,i+2*j-2)
      plot(devals{i},ddevals{j}-fdevals{j})
   end
end
    

%% Perform the Profiled Estimation

[newpars,newDEfd] = Profile_GausNewt(startpars,lsopts_out,tDEfd,fn,...
    lambda,Ycell,Tcell,wts,[],lsopts_in);

disp(['New parameter estimates: ',num2str(newpars')]);

% plot smooth with profile-estimated parameters

figure(6)
devals = eval_fdcell(tfine,newDEfd,0);
for i = 1:length(Ycell)
    subplot(length(Ycell),1,i)
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    hold off
    if i==1
        ylabel('\fontsize{13} V','rotation',0)
%        title(['\fontsize{13} Raw data (.), ', ...
%               'profiled solution (r-) and true path (g-)'])
    else
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} R','rotation',0)
    end
end


% Now try diagnostics


fn.fn       = @forcingfun;        % RHS function

fn.dfdx     = @forcingdfdx;       % Derivative wrt inputs (Jacobian)
fn.dfdp     = @forcingdfdp;       % Derviative wrt parameters

fn.d2fdx2   = @forcingd2fdx2;     % Hessian wrt inputs
fn.d2fdxdp  = @forcingd2fdxdp;    % Cross derivatives wrt inputs and parameters
fn.d2fdp2   = @forcingd2fdp2;     % Hessian wrt parameters

fn.d3fdx2dp = @forcingd3fdx2dp;   % Third derivative wrt inputs, inputs, pars 
fn.d3fdx3   = @forcingd3fdx3;     % Third derivative wrt inputs
fn.d3fdxdp2 = @forcingd3fdxdp2;   % Third derivative wrt inputs, pars and pars

% Define the original equations

fn_extras.fn      = @runzlefun;    % Original function
fn_extras.dfdx    = @runzledfdx;   % First derivative
fn_extras.d2fdx2  = @runzled2fdx2; % Second derivative

fn_extras.extras  = [];        % Original information to fn_extras.fn. 

% We also need to add parameters:

fn_extras.pars = newpars;

pnknots = 251;
pnbasis = pnknots + norder - 2;        
pknots = linspace(range(1),range(2),pnknots);
pbasis_obj = MakeBasis(range,pnbasis,norder,pknots,quadvals,3);

fn_extras.basisp = {pbasis_obj};   % Cell array of basis functions.
fn_extras.which  = 1;          


pen   = @forcingpen;   % Penalty
dpen  = @forcingdpen;  % Derivative with respect to forcing co-efficients
d2pen = @forcingd2pen; % Second derivative

pen_extras. basis = fn_extras.basisp; % Same basis as fn_extras. 
pen_extras.deg    = [0 2];           % Penalize the second derivative
pen_extras.lambda = [1e-4 1e-4];      % Smoothing parameter.                            

startpars = zeros(getnbasis(pbasis_obj),1);

% Now do the optimization and hope it doesn't fall over

[fcoefs,smooths] = Profile_GausNewt(startpars,lsopts_out,newDEfd,fn,lambda,...
    Ycell,Tcell,[],[],lsopts_in,fn_extras,pen,dpen,pen_extras);                                 

forces = Make_fdcell(fcoefs,{pbasis_obj});


% Evaluate these pointwise so their values can be compared. 

t = 0:1000;

ss = eval_fdcell(t,smooths);
fs = eval_fdcell(t,forces);

% First of all plot the smooth with the data; we observe a reasonable
% correspondence. 

figure(3)
for i = 1:2
    subplot(2,1,i);
    plot(t,ss{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    hold off
end


figure(4)
for i = 1:length(fs)
    subplot(length(fs),1,i)
    plot(t,fs{i})
end

%% Diagnostics
%
% To try to evaluate how we should change the system in order to provide a
% more accurate fit, we plot each of the forcing functions against each of
% the components and observe whether there appears to be a systematic
% relationship. 

figure(5)
k = 0;
for i = 1:length(ss)
    for j = 1:length(fs)
        k = k+1;
        subplot(length(ss),length(fs),k)
        plot(ss{i},fs{j})
    end
end

figure(6)
plot3(ss{1},ss{2},fs{1})



% Model-Building for spiking rate of a Zebra Fish/

% RHS Functions

fn.fn       = @genlinfun;       % RHS function

% Now derivatives of the function with respect to system components and
% parameters 

fn.dfdx     = @genlindfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @genlindfdp;      % Derviative wrt parameters

% Now we need functions to compute all three sets of second derivatives:

fn.d2fdx2   = @genlind2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @genlind2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @genlind2fdp2;    % Hessian wrt parameters.    


                          
startpars = zeros(4,1);

%% Load the data

load 'C:\Documents and Settings\Giles\Desktop\zebrafish\test03.asc'

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

lambda = 1e4;      
lambda = lambda * wts;

lambda0 = 1;        

nknots = 501;    
norder = 6;
nquad = 5;     

%% Profiling optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-8,'TolX',1e-10);

lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


%% Setting up Functional Data Objects

range = [0 1000];

const_basis = create_constant_basis([0 1000]);
const_fd = fd(1,const_basis);
fn_extras.force =  {const_fd};

fn_extras.force_sub = [2 1];

fn_extras.sub = [2 1; 2 2];
fn_extras.mat = zeros(2,2);
fn_extras.mat(1,2) = 1;


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
plot(Tcell{i},Ycell{1},'b.');
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

%% Perform the Profiled Estimation

[newpars,newDEfd] = Profile_GausNewt(startpars,lsopts_out,DEfd1,fn,...
    lambda,Ycell,Tcell,wts,[],lsopts_in,fn_extras);

disp(['New parameter estimates: ',num2str(newpars')]);

% plot smooth with profile-estimated parameters

basis_cell = getcellbasis(DEfd1);
[newcoefs,newDEfd] = genlin_smooth(Ycell,Tcell,wts,basis_cell,lambda,...
    newpars,[],fn_extras);

figure(4)
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


% Now do the linear approximation

A = fn_extras.mat;
A(2,1:2) = newpars(1:2);
const_force = {fd(0,const_basis); fd(newpars(3),const_basis)};

[smooths,forces] = linforceest(basis_cell,basis_cell(1),A,1,10000,...
    1,2,Tcell,Ycell,wts,0,@eval_fdcell2,const_force);


ss = eval_fdcell(tfine,smooths);
fs = eval_fdcell(tfine,forces);

% First of all plot the smooth with the data; we observe a reasonable
% correspondence. 

figure(3)
axes('fontsize',14)
subplot(2,1,1)
plot(tfine,Ycell{1},'.');
hold on
plot(tfine,ss{1},'r','linewidth',2);
hold off

subplot(2,1,2)
plot(tfine,ss{2},'r','linewidth',2);


figure(4)
plot(tfine,fs{1},'linewidth',2);

% Now do the diagnostics

figure(5)
axes('fontsize',14)

subplot(2,1,1)
plot(ss{1},fs{1},'linewidth',2)
subplot(2,1,2)
plot(ss{2},fs{1},'linewidth',2)


figure(6)
plot3(ss{1},ss{2},fs{1})


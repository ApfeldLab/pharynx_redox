% Model-Building for spiking rate of a Zebra Fish/

% RHS Functions

fn.fn       = @fhnfun;       % RHS function

% Now derivatives of the function with respect to system components and
% parameters 

fn.dfdx     = @fhndfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @fhndfdp;      % Derviative wrt parameters

% Now we need functions to compute all three sets of second derivatives:

fn.d2fdx2   = @fhnd2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @fhnd2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @fhnd2fdp2;    % Hessian wrt parameters.    


startpars = [0 0 1]';

%% Load the data

load 'zebrafish\test03.asc'

Ycell = cell(1,2);
Tcell = cell(1,2);

Ycell{1} = test03(2000:3000,1);
Tcell{1} = 0.015*(0:1000);

tfine = Tcell{1};

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
    'Display','iter','MaxIter',100,'TolFun',1e-8,'TolX',1e-10);

lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',100,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',100,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


%% Setting up Functional Data Objects

range = [tfine(1) tfine(1001)];

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



%% Perform the Profiled Estimation

[newpars,newDEfd] = Profile_GausNewt(startpars,lsopts_out,DEfd1,fn,...
    lambda,Ycell,Tcell,wts,[],lsopts_in);

disp(['New parameter estimates: ',num2str(newpars')]);

% plot smooth with profile-estimated parameters

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






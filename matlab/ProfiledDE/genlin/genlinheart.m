%% Profile Estimation Experiments - General Forced Linear Systems

% Load and smooth the data

load('../heartprint.txt');

t = heartprint(1:1000,1);
V = heartprint(1:1000,2);
R = heartprint(1:1000,3);

Ycell = {V R};

knots = linspace(t(1),t(1000),250);
norder = 6;
nquad = 5;

quadpoints = MakeQuadPoints(knots,nquad);
basis_obj = MakeBasis([t(1) t(1000)],length(knots)+norder-2,norder,...
    knots,quadpoints,4);
Lfd_obj = fdPar(basis_obj,4,1e-15);

Lfd_cell = {Lfd_obj Lfd_obj};

smooth = smoothfd_cell(t,Ycell,Lfd_cell};
ss = eval_fdcell(t,smooth);

figure(1)
plot(t,V,'b')
hold on
plot(t,R,'g')
plot(t,ss{1},'r')
plot(t,ss{2},'r')
title('\fontsize{20} HeartPrint and Smooth')
xlabel('\fontsize{20} seconds')

% Set up a 2nd-order system

odefn       = @genlinfunode;    % Function for ODE solver (exact)


fn.fn       = @genlinfun;       % RHS function
fn.dfdx     = @genlindfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @genlindfdp;      % Derviative wrt parameters

fn.d2fdx2   = @genlind2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @genlind2fdxdp;   % Hessian wrt inputs and parameters
fn.d2fdp2   = @genlind2fdp2;    % Hessian wrt parameters.    

fn.d3fdx3   = @genlind3fdx3;    % Third derivative wrt inputs.
fn.d3fdx2dp = @genlind3fdx2dp;  % Third derivative wrt intputs, inputs and pars.
fn.d3fdxdp2 = @genlind3fdxdp2;  % Third derivative wrt inputs, pars and pars. 
                                % dimensions = time, component, input,
                                % parameters                     
                          


more = [];
more.mat = zeros(4,4);
more.mat(1:2,1:2) = eye(2);

more.sub = [3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];         


Tcell = {t t [] []};
Ycell = {V R [] []};

basis_cell = cell(1,4);
basis_cell(:) = {basis_obj};

Lfd_cell = cell(4,1);
Lfd_cell(:) = {Lfd_obj};

%% Very rough and ready parameters

dss = cell2mat(eval_fdcell(t,smooth,1));
d2ss = cell2mat(eval_fdcell(t,smooth,2));

X = [cell2mat(ss) dss];

startpars = (X'*X)\X'*d2ss;

disp('Starting Parameter Estimates')
disp(startpars)

plot(d2ss - X*startpars)

startpars = [startpars(:,1); startpars(:,2)];


%% Fitting parameters

lambda  = 1e4;   % Smoothing for model-based penalty

%% Profiling optimisation control

maxit1 = 1000;      % Maximum iterations interior of profiling
maxit0 = 1000;      % Maximum iterations for outer optimization

lsopts_out = optimset('DerivativeCheck','on','Jacobian','on',...
    'Display','iter','MaxIter',maxit0,'TolFun',1e-8,'TolX',1e-10);

% Other observed optimiation control
lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',maxit0,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

% Optimiation control within profiling
lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',maxit1,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);



%% Smooth the data

DEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);
coefs = getcellcoefs(DEfd);

devals = eval_fdcell(t,DEfd,0);
for i = 1:length(Ycell)
    subplot(length(Ycell),1,i);
    plot(t,devals{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    hold off;
end


% Call the Gauss-Newton solver

[newcoefs,tDEfd] = genlin_smooth(Ycell,Tcell,[],basis_cell,lambda,startpars,[],more);

% Plot results along with exact solution

devals = eval_fdcell(t,tDEfd,0);
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(t,devals{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    hold off
end



%% Perform the Profiled Estimation

[newpars,newDEfd_cell] = Profile_GausNewt(startpars,lsopts_out,tDEfd,fn,...
    lambda,Ycell,Tcell,[],[],lsopts_in,more);

disp('New parameter values:')
disp(newpars');



%% Plot Smooth with Profile-Estimated Parameters

devals = eval_fdcell(tfine,newDEfd_cell,0);
for i = 1:length(path)
    subplot(length(path),1,i)
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
    plot(Tcell{i},Ycell{i},'b.');
    plot(plot_time,plot_path(:,i),'c');
    hold off
end



%% Comparison with Smooth Using True Parameters

[truecoefs,trueDEfd_cell] = genlin_smooth(path,Tcell,[],basis_cell,...
    lambda,pars,[],more);

devals = eval_fdcell(tfine,trueDEfd_cell,0);
for i = 1:length(path)
    subplot(length(path),1,i)
    plot(plot_time,plot_path(:,i),'c')
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
    plot(plot_time,plot_path(:,i),'c');
    plot(Tcell{i},Ycell{i},'b.');
    hold off;
end



%% Squared Error Performance

% Squared error for estimated parameters

newpreds = eval_fdcell(Tcell,newDEfd_cell,0);
new_err = cell(length(newpreds));
for i = 1:length(path)
    new_err{i} = wts(i)*(newpreds{i} - Ycell{i}).^2;
end

new_err = mean(cell2mat(new_err));

% Squared error for true parameters


truepreds = eval_fdcell(Tcell,trueDEfd_cell,0);
true_err = cell(length(truepreds));
for i = 1:length(path)
    true_err{i} = wts(i)*(truepreds{i} - Ycell{i}).^2;
end

true_err = mean(cell2mat(true_err));

% print out a comparision

disp(['Estimated sqrd error: ',num2str(new_err)])
disp(['True sqrd error:      ',num2str(true_err)]);


%% Calculate Sample Information and Variance-Covariance Matrices

% Hessian of squared error with respect to parameters

d2Jdp2 = make_d2jdp2(newDEfd_cell,fn,Ycell,Tcell,lambda,newpars,[],wts,more);

% Second derivatives with respect to parameters and observations

d2JdpdY = make_d2jdpdy(newDEfd_cell,fn,Ycell,Tcell,lambda,newpars,[],...
    wts,more);

% Resulting derivative of parameters with respect to observations

dpdY = -d2Jdp2\d2JdpdY;

% Variance of observations:

S = make_sigma(DEfd,Tcell,Ycell,0);

% Resulting parameter covariance matrix:

Cov = dpdY*S*dpdY';

%  Standard errors

StdDev = sqrt(diag(Cov));

%  Correlations

Corr = Cov./(StdDev*StdDev');

%  Display these results

disp('Approximate covariance matrix for parameters:')
disp(num2str(Cov))

disp('Approximate standard errors of parameters:')
disp(num2str(StdDev'))

disp('Approximate correlation matrix for parameters:')
disp(num2str(Corr))


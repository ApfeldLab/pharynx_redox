%% Simulation of Parameter Estimation for Perturbed Linear Systems

odefn    = @genlinfunode;    % Function for ODE solver (exact)

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
                                                   
%% Various Parameters

more = [];

y0 = [-1,1];                      
pars = [-1; 2; -1; 1]; 

sigma = 1;                    

jitter = 0.2;                   
nrep = 600;

%% Observation times

tspan = 0:0.05:20;  
obs_pts{1} = 1:length(tspan);       
obs_pts{2} = 1:length(tspan);      


%% Fitting parameters


lambdas = 10.^(-6:0.5:2);     
lambda0 = 1;        


nknots = 401;
norder = 6;
nquad = 5;     

%% Profiling optimisation control

lsopts_out = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','iter','MaxIter',1000,'TolFun',1e-8,'TolX',1e-10);

lsopts_in = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','off','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);

lsopts_other = optimset('DerivativeCheck','off','Jacobian','on',...
    'Display','on','MaxIter',1000,'TolFun',1e-14,'TolX',1e-14,...
    'JacobMult',@SparseJMfun);


%% Setting up Functional Data Objects

range = [min(tspan),max(tspan)];

knots_cell = cell(1,length(y0));
knots_cell(:) = {linspace(range(1),range(2),401)};

basis_cell = cell(1,length(y0)); 
Lfd_cell = cell(1,length(y0));

nbasis = zeros(length(y0),1);

bigknots = knots_cell{1};               
nbasis(1) = length(knots_cell{1}) + norder - 2;          

for(i = 2:length(y0))
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(bigknots,nquad);   

for(i = 1:length(y0))
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  
        knots_cell{i},quadvals,4);                        
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         
end

%% Calculate paths

pbasis = create_bspline_basis([0,20],3,0,[0,7,14,20]);

fd_obj = fd([0; 2; 0],pbasis);
more.force = {fd_obj};
more.force_mat = [1; 0];

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,[pars; 1; 0],more);


figure(1)
plot(full_time,full_path,'LineWidth',2)

figure(2)
plot(full_path(:,1),full_path(:,2))

%% Set up observations

Tcell = cell(1,size(full_path,2));
path = Tcell;

for(i = 1:length(obs_pts))
    Tcell{i} = full_time(obs_pts{i});
    path{i} = full_path(obs_pts{i},i);
end

wts = [];  
for(i = 1:length(path))
    wts(i) = 1./sqrt(var(path{i}));
end


Ycell = path;
for(i = 1:length(path))
    Ycell{i} = path{i} + sigma*randn(size(path{i}));
end

% Smooth the data

sse = Inf;
for(i = 1:length(lambdas))
    for(j = 1:length(Ycell))
        Lfd_cell{j} = fdPar(basis_cell{j},1,lambdas(i));
    end
    
    tDEfd = smoothfd_cell(Ycell,Tcell,Lfd_cell);

    tsse = sum( sum( (cell2mat(path)-cell2mat(eval_fdcell(Tcell,tDEfd))).^2 ));
    
    if tsse < sse
        sse = tsse;
        DEfd = tDEfd;
        disp(i)
    else if i>1
            break;
        end
    end
end    
    
A = reshape(pars,2,2)';

sse = Inf;
for(i = 1:length(lambdas))

    [tsmooths,tforces] = linforceest(basis_cell,basis_cell(1),A,1,10000,...
        lambdas(i),1,Tcell,Ycell,wts);

    tss = cell2mat(eval_fdcell(Tcell,tsmooths));
    tsse = sum(sum( (tss-cell2mat(path)).^2 ));
    
    if(tsse<sse)
        sse = tsse;
        smooths = tsmooths;
        forces = tforces;
        disp(i)
    else
        if i>1
            break;
        end
    end
end
        
% Evaluate these pointwise so their values can be compared. 

tfine = full_time;
ss = eval_fdcell(tfine,smooths);
fs = eval_fdcell(tfine,forces);

% First of all plot the smooth with the data; we observe a reasonable
% correspondence. 

figure(3)
for i = 1:length(path)
    subplot(length(path),1,i);
    plot(Tcell{i},Ycell{i},'b.');
    hold on;    
    plot(full_time,full_path(:,i),'c','LineWidth',2);
    plot(tfine,ss{i},'r--','LineWidth',2);
    hold off
    if i==1
        ylabel('\fontsize{20} X')
%        title(['\fontsize{13} Raw data (.), ', ...
%               'smoothing solution (r-) and true path (g-)'])
    else
        xlabel('\fontsize{20} t')
        ylabel('\fontsize{20} Y')
    end
end

% Now plot estimated forcing functions

fvals = eval_fd(tfine,fd_obj);

vals = cell2mat(eval_fdcell(tfine,DEfd));
svals = cell2mat(eval_fdcell(tfine,DEfd,1));

disc = svals - vals*A';

figure(4)
plot(tfine,fvals,'LineWidth',2)
hold on
plot(tfine,fs{1},'r--','LineWidth',2)
plot(tfine,disc(:,1),'g-.','LineWidth',2)
hold off
xlabel('\fontsize{20} t')
ylabel('\fontsize{20} g')




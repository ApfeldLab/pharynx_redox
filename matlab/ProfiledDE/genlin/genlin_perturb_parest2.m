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

sigma = 0.5;                    

jitter = 0.2;                   
nrep = 600;

%% Observation times

tspan = 0:0.05:20;  
obs_pts{1} = 1:length(tspan);       
obs_pts{2} = 1:length(tspan);      


%% Fitting parameters


lambdas = 10.^(-2:0.2:4);     
lambda0 = 1;        


nknots = 401;
norder = 3;
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
        knots_cell{i},quadvals,1);                        
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         
end

%% Calculate paths

pcoef = zeros(nbasis(1),1);
pcoef(floor(nbasis(1)/2):(ceil(nbasis(1)/2)+5)) = 10;
fd_obj = fd(pcoef,basis_cell{1});
more.force = {fd_obj};
more.force_mat = [1; 0];

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,[pars; 1; 0],more);

plot(full_time,full_path,'LineWidth',2)

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


% Smooth the data

DEfd = smoothfd_cell(path,Tcell,Lfd_cell);
coefs = getcellcoefs(DEfd);


% Start the simulations

errs = zeros([nrep length(lambdas) 3]);
pppars = zeros(nrep,length(pars));

Ycell = path;
for(i = 1:length(path))
    Ycell{i} = path{i} + sigma*randn(size(path{i}));
end

for(g = 1:nrep)

    ppars = pars + jitter*randn(length(pars),1);
    
    if(g == 1)
        ppars = pars;
    end
    
    disp(ppars')    
    pppars(g,:) = ppars;
    
    for(h = 1:length(lambdas))
        disp([g h])

        lambda = lambdas(h) * wts;
        disp(lambda)

        %% Comparison with Smooth Using True Parameters

        coefs = getcellcoefs(DEfd);

        [coefs,DEfd_cell] = genlin_smooth(Ycell,Tcell,wts,basis_cell,...
            lambda,ppars,[],[]);
        

        newpreds = eval_fdcell(Tcell,DEfd_cell,0);

        new_err = cell(length(newpreds));
        new_err_true = cell(length(newpreds));
        for(i = 1:length(path))
            new_err{i} = wts(i)*(newpreds{i} - Ycell{i}).^2;
            new_err_true{i} = wts(i)*(newpreds{i} - path{i}).^2;
        end

        new_err = mean(cell2mat(new_err));
        new_err_true = mean(cell2mat(new_err_true));

        % Evaluate size of penalty
            
        dpreds = cell2mat(eval_fdcell(quadvals(:,1),DEfd_cell,1));
        preds = cell2mat(eval_fdcell(quadvals(:,1),DEfd_cell,0));
        
        pparsM = reshape(ppars,2,2);
        
        pred_err = sum(quadvals(:,2).*sum(dpreds-preds*pparsM',2).^2);
        
        % Record
        
        disp([g h new_err new_err_true pred_err])
        
        errs(g,h,:) = [new_err new_err_true pred_err];
       
    end

end

save 'genlin_perturb_parest3.mat' errs pppars

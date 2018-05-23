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


%% Observation times

tspan = 0:0.05:20;  
obs_pts{1} = 1:length(tspan);       
obs_pts{2} = 1:length(tspan);      


%% Fitting parameters


lambdas = 10.^8;
lambdap = 10.^(-7:2:-1);     
levels = 0:0.1:0.5;

nrep = 10000;
nknots = 51;
norder = 6;
nquad = 5;     


%% Setting up Functional Data Objects

range = [min(tspan),max(tspan)];

knots_cell = cell(1,length(y0));
knots_cell(:) = {linspace(range(1),range(2),401)};

basis_cell = cell(1,length(y0)); 
Lfd_cell = cell(1,length(y0));

nbasis = zeros(length(y0),1);

bigknots = knots_cell{1};               
nbasis(1) = length(knots_cell{1}) + norder - 2;          

for i = 2:length(y0)
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(bigknots,nquad);   

for i = 1:length(y0)
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  
        knots_cell{i},quadvals,4);                                
end

A = reshape(pars,2,2)';

pbasis = create_bspline_basis([0,20],3,0,[0,7,14,20]);
fd_obj = fd([0; 1; 0],pbasis);

wts = [1 1];

%% Set up iterations

wald = zeros(length(levels),length(lambdas),length(lambdap),nrep);
lr = zeros(length(levels),length(lambdas),length(lambdap),nrep);

sse1 = zeros(length(levels),length(lambdas),nrep);
df1 = zeros(length(lambdas),1);

Slam1 = zeros(sum(nbasis),2*length(tspan),length(lambdas));

sse2 = zeros(length(levels),length(lambdas),length(lambdap),nrep);
df2 = zeros(length(lambdas),length(lambdap));

fCov = zeros(nbasis(1),2*length(tspan),length(lambdas),length(lambdap));
fCov2 = zeros(nbasis(1),nbasis(1),length(lambdas),length(lambdap));

cCov = zeros(sum(nbasis),2*length(tspan),length(lambdas),length(lambdap));

fcoef = zeros(nbasis(1),length(lambdas),length(lambdap),length(levels),nrep);

for h = 1:length(levels)

    more.force = {fd_obj};
    more.force_mat = [1; 0];

    odeopts = odeset('RelTol',1e-13);
    [full_time,full_path] = ode45(odefn,tspan,y0,odeopts,[pars; levels(h); 0],more);
    
    for g = 1:nrep


        Tcell = cell(1,size(full_path,2));
        path = Tcell;

        for i = 1:length(obs_pts)
            Tcell{i} = full_time(obs_pts{i});
            path{i} = full_path(obs_pts{i},i);
        end


        Ycell = path;
        for  i = 1:length(path)
            Ycell{i} = path{i} + sigma*randn(size(path{i}));
        end

        y = cell2mat(Ycell');
        
        % Smooth the data

        for i = 1:length(lambdas)
            
            if Slam1(1,1,i) == 0
            
               [coefs,DEfd,Slam1(:,:,i),df1(i)] = genlin_smooth(Ycell,Tcell,wts,basis_cell,...
                    lambdas(i),pars);

            else
                
                coefs = Slam1(:,:,i)*y;
                DEfd = Make_fdcell(coefs,basis_cell);
                
            end
            
            sse1(h,i,g)=sum(sum((y-cell2mat(eval_fdcell(Tcell,DEfd)')).^2));

            for j = 1:length(lambdap)
                
                disp([h g i j])
                
                if cCov(1,1,i,j) == 0
                
                    [smooths,forces,fCov(:,:,i,j),cCov(:,:,i,j),df2(i,j)] = linforceest(basis_cell,...
                        basis_cell(1),A,1,lambdas(i),lambdap(j),1,Tcell,Ycell,wts);
                    
                    fCov2(:,:,i,j) = inv(fCov(:,:,i,j)*fCov(:,:,i,j)');
        
                    tfcoef = getcellcoefs(forces);
                else
                   
                   tfcoef = fCov(:,:,i,j)*y;
                   scoef = cCov(:,:,i,j)*y;
                   smooths = Make_fdcell(scoef,basis_cell);
                    
                end
                    
                
                fcoef(:,i,j,h,g) = tfcoef;
                S = make_sigma(smooths,Tcell,Ycell,0);

                sigmahat(i,j,h,g) = S(1,1);

                wald(h,i,j,g) = tfcoef'*fCov2(:,:,i,j)*tfcoef;

                sse2(h,i,j,g)=sum(sum((y-cell2mat(eval_fdcell(Tcell,smooths)')).^2));

                lr(h,i,j,g) = sse1(h,i,g) - (sse2(h,i,j,g));

            end
        end
    end
    save 'linforce_noparest_thisdate.mat'
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
        pwald(:,i,j) = sum(squeeze(wald(:,i,j,:))>wald05(j,i),2);
        plr(:,i,j) = sum(squeeze(lr(:,i,j,:))>lr05(j,i),2);
    end
end
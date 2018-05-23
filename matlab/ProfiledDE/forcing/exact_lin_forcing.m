
% Basis data setup

pars = [-1; 2; -1; 1]; 
A = reshape(pars,2,2)';
y0 = [-1 1];

tspan = 0:0.05:20;  
obs{1} = 1:length(tspan);       
obs{2} = 1:length(tspan);    

Tcell = {tspan(obs{1}), tspan(obs{2})};

sigma = 1;        

% ODE Solver

pbasis = create_bspline_basis([0,20],3,0,[0,7,14,20]);
fd_obj = fd([0; 1; 0],pbasis);

addpath('../genlin')
odefn    = @genlinfunode;  
odeopts = odeset('RelTol',1e-13);

% Experiment parameters
           
lambdap = 10.^(-7:2:3);     
levels = 0:0.1:0.5;

nrep = 10000;

SSE1 = zeros(length(levels),nrep);

SSEa = zeros(length(levels),nrep,length(lambdap));
SSEa1 = zeros(length(levels),nrep,length(lambdap));
SSEa2 = zeros(length(levels),nrep,length(lambdap));

SSEs = zeros(length(levels),nrep,length(lambdap));
SSEs1 = zeros(length(levels),nrep,length(lambdap));
SSEs2 = zeros(length(levels),nrep,length(lambdap));

lra = zeros(length(levels),nrep,length(lambdap));
lra1 = zeros(length(levels),nrep,length(lambdap));
lra2 = zeros(length(levels),nrep,length(lambdap));

lrs = zeros(length(levels),nrep,length(lambdap));
lrs1 = zeros(length(levels),nrep,length(lambdap));
lrs2 = zeros(length(levels),nrep,length(lambdap));

% Define the bases

nknots = 51;
norder = 6;
nquad = 5;   

range = [min(tspan),max(tspan)];

knots_cell = cell(1,length(y0));
knots_cell(:) = {linspace(range(1),range(2),nknots)};

basis_cell = cell(1,length(y0)); 
Lfd_cell = cell(1,length(y0));

nbasis = zeros(length(y0),1);

bigknots = knots_cell{1};               
nbasis(1) = length(knots_cell{1}) + norder - 2;          

for i = 2:length(y0)
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(tspan,nquad);   
quadpts = unique(quadvals(:,1));
delta = quadpts(2)-quadpts(1);
quadvals = [quadpts delta*(ones(size(quadpts)))];

for i = 1:length(y0)
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  
        knots_cell{i},quadvals,4);                                
end

% And a penalty matrix

D = cell(3,1);
D{1} = zeros(2);
D{2} = eval_penalty(basis_cell{1},1);
D{3} = eval_penalty(basis_cell{2},1);

Da = mattdiag_cell(D,0);
D1 = mattdiag_cell(D(1:2),0);
D2 = mattdiag_cell(D([1 3]),0);

% Now the terms in the linear expansion

T = tspan;
qT = 1:(nquad-1):size(quadvals(:,1));

for i = 1:length(T)
    E(:,:,i) = expm(A*T(i));
end

X = [squeeze(E(1,:,obs{1}))'; squeeze(E(2,:,obs{2}))'];

for i = 1:length(quadpts)
   Q(:,:,i) = expm(-A*quadpts(i));
end

Zmat = getvalues_cell(basis_cell,0);

for i = 1:2
       nQ(i,:,:) = cumsum( [diag(squeeze(Q(i,1,:)))*Zmat{1} ...
                           diag(squeeze(Q(i,2,:)))*Zmat{2}] )*delta;
end

nnQ1(1,:,:) = diag(squeeze(E(1,1,:)))*squeeze(nQ(1,qT,:)); 
nnQ2(1,:,:) = diag(squeeze(E(1,2,:)))*squeeze(nQ(2,qT,:));
nnQ1(2,:,:) = diag(squeeze(E(2,1,:)))*squeeze(nQ(1,qT,:));
nnQ2(2,:,:) = diag(squeeze(E(2,2,:)))*squeeze(nQ(2,qT,:));

nnQ = nnQ1+nnQ2;

Z = [squeeze(nnQ(1,obs{1},:)); squeeze(nnQ(2,obs{2},:))];

Zs = mattdiag_cell(eval_basis_cell(Tcell,basis_cell,0),0);

XX = [X Z];
XXs = [X Zs];

XX1 = [X Z(:,1:nbasis(1))];
XX2 = [X Z(:,nbasis(1)+(1:nbasis(2)))];

XXs1 = [X Zs(:,1:nbasis(1))];
XXs2 = [X Zs(:,nbasis(1)+(1:nbasis(2)))];

% Finally, the smoother matrices

S1 = X*inv(X'*X)*X';

I = eye(length(obs{1}) + length(obs{2}));
df1 = trace((I-S1)'*(I-S1));

for i = 1:length(lambdap)

    Sa(:,:,i) = XX*inv(XX'*XX + lambdap(i)*Da)*XX';
    Sa1(:,:,i) = XX1*inv(XX1'*XX1 + lambdap(i)*D1)*XX1';
    Sa2(:,:,i) = XX2*inv(XX2'*XX2 + lambdap(i)*D2)*XX2';
    
    Ss(:,:,i) = full(XXs*inv(XXs'*XXs + lambdap(i)*Da)*XXs');
    Ss1(:,:,i) = full(XXs1*inv(XXs1'*XXs1 + lambdap(i)*D1)*XXs1');
    Ss2(:,:,i) = full(XXs2*inv(XXs2'*XXs2 + lambdap(i)*D2)*XXs2');
    
    
    dfa(i) = trace( (I-Sa(:,:,i))'*(I-Sa(:,:,i)) );
    dfa1(i) = trace( (I-Sa1(:,:,i))'*(I-Sa1(:,:,i)) );
    dfa2(i) = trace( (I-Sa2(:,:,i))'*(I-Sa2(:,:,i)) );

    dfs(i) = trace( (I-Ss(:,:,i))'*(I-Ss(:,:,i)) );
    dfs1(i) = trace( (I-Ss1(:,:,i))'*(I-Ss1(:,:,i)) );
    dfs2(i) = trace( (I-Ss2(:,:,i))'*(I-Ss2(:,:,i)) );

    qfa(i) = chisqq(0.95,df1-dfa(i));
    qfa1(i) = chisqq(0.95,df1-dfa1(i));
    qfa2(i) = chisqq(0.95,df1-dfa2(i));

    qfs(i) = chisqq(0.95,df1-dfs(i));
    qfs1(i) = chisqq(0.95,df1-dfs1(i));
    qfs2(i) = chisqq(0.95,df1-dfs2(i));
end

for i = 1:length(lambdap)
    dfa_2(i) = trace(Sa(:,:,i));
    dfa1_2(i) = trace(Sa1(:,:,i));
    dfa2_2(i) = trace(Sa2(:,:,i));
    
    dfs_2(i) = trace(Ss(:,:,i));
    dfs1_2(i) = trace(Ss1(:,:,i));
    dfs2_2(i) = trace(Ss2(:,:,i));


    qfa_2(i) = chisqq(0.95,dfa_2(i)-2);
    qfa1_2(i) = chisqq(0.95,dfa1_2(i)-2);
    qfa2_2(i) = chisqq(0.95,dfa2_2(i)-2);

    qfs_2(i) = chisqq(0.95,dfs_2(i)-2);
    qfs1_2(i) = chisqq(0.95,dfs1_2(i)-2);
    qfs2_2(i) = chisqq(0.95,dfs2_2(i)-2);

end




% Now run the simulation

for i = 1:length(levels)
    
    more.force = {fd_obj};
    more.force_mat = [1; 0];

    [full_time,full_path] = ode45(odefn,tspan,y0,odeopts,[pars; levels(i); 0],more);

    fY = [full_path(obs{1},1); full_path(obs{2},2)];
    
    for j = 1:nrep
        
        Y = fY + sigma*randn(size(fY));
        SSE1(i,j) = sum( (Y - S1*Y).^2 );
        
        for k = 1:length(lambdap)
            disp([i j k])

           SSEa(i,j,k) = sum( (Y - Sa(:,:,k)*Y).^2 );
           SSEa1(i,j,k) = sum( (Y - Sa1(:,:,k)*Y).^2 );
           SSEa2(i,j,k) = sum( (Y - Sa2(:,:,k)*Y).^2 );
           
           SSEs(i,j,k) = sum( (Y - Ss(:,:,k)*Y).^2 );
           SSEs1(i,j,k) = sum( (Y - Ss1(:,:,k)*Y).^2 );
           SSEs2(i,j,k) = sum( (Y - Ss2(:,:,k)*Y).^2 );
           
           lra(i,j,k) = SSEa(i,j,k) - SSE1(i,j);
           lra1(i,j,k) = SSEa1(i,j,k) - SSE1(i,j);
           lra2(i,j,k) = SSEa2(i,j,k) - SSE1(i,j);
           
           lrs(i,j,k) = SSEs(i,j,k) - SSE1(i,j);
           lrs1(i,j,k) = SSEs1(i,j,k) - SSE1(i,j);
           lrs2(i,j,k) = SSEs2(i,j,k) - SSE1(i,j);
       end        
        
    end    
end


for i = 1:length(lambdap)
    
    plra(:,i) = mean( -lra(:,:,i) > qfa_2(i),2 );
    plra1(:,i) = mean( -lra1(:,:,i) > qfa1_2(i),2 );
    plra2(:,i) = mean( -lra2(:,:,i) > qfa2_2(i),2 );

    plrs(:,i) = mean( -lrs(:,:,i) > qfs_2(i),2 );
    plrs1(:,i) = mean( -lrs1(:,:,i) > qfs1_2(i),2 );
    plrs2(:,i) = mean( -lrs2(:,:,i) > qfs2_2(i),2 );
end




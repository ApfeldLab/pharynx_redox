% Try just using smoothing methods:


%% Load the data

load 'zebrafish\test03.asc'

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


%% Fitting parameters

lambda0 = 1;        

nknots = 501;    
norder = 6;
nquad = 5;     


%% Setting up Functional Data Objects

range = [0 1000];

nbasis = nknots + norder - 2;        
knots = linspace(range(1),range(2),nknots);

quadvals = MakeQuadPoints(knots,nquad);   

basis_obj = MakeBasis(range,nbasis,norder,knots,quadvals,3);
Lfd_obj = fdPar(basis_obj,3,lambda0);

Lfd_cell = {Lfd_obj Lfd_obj};

basis_cell = {basis_obj basis_obj};

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


% We also want to estimate parameters from this smooth

ddevals = eval_fdcell(tfine,DEfd1,1);

figure(3)
axes('fontsize',14)
plot(devals{1},ddevals{1},'linewidth',2);
xlabel('\fontsize{20} V')
ylabel('\fontsize{20} dV/dt','rotation',90)

figure(4)
d2devals = eval_fdcell(tfine,DEfd1,2);
plot3(devals{1},ddevals{1},d2devals{1});
xlabel('\fontsize{20} V')
ylabel('\fontsize{20} dV')
zlabel('\fontsize{20} d2V')


% Now various smooths

x1 = devals{1};
x2 = ddevals{1};
y = d2devals{1};
o1 = ones(size(x1));

X = [o1 x1 x2];

par1 = (X'*X)\X'*y;

res1 = y - X*par1;

plot3(x1,x2,res1);


% Add an interaction

X = [X x1.*x2];

par2 = (X'*X)\X'*y;

res2 = y - X*par2;

plot3(x1,x2,res2);
xlabel('x1')
ylabel('dx1')

% Looks vaguely cubic:

X = [X x1.^2 x1.^3];

par3 = (X'*X)\X'*y;

res3 = y - X*par3;

plot3(x1,x2,res3);
xlabel('x1')
ylabel('dx1')

% Interactions with nonlinear terms

X = [X x2.*x1.^2 x2.*x1.^3];

par4 = (X'*X)\X'*y;

res4 = y - X*par4;

plot3(x1,x2,res4);
xlabel('x1')
ylabel('dx1')


% Try adding nonlinearities in derivatives:

X = [X x2.^2.*x1 x2.^2.*x1.^2];

par5 = (X'*X)\X'*y;

res5 = y - X*par5;

plot3(x1,x2,res5);
xlabel('x1')
ylabel('dx1')

% Ok, full fourth order

X = [X x2.^3.*x1 x2.^2 x2.^3  x2.^4 x1.^4];

par6 = (X'*X)\X'*y;

res6 = y - X*par6;

axes('fontsize',14)
plot3(x1,x2,res6);
xlabel('\fontsize{20} v')
ylabel('\fontsize{20} dv')
zlabel('\fontsize{20} d2v')

% Lets look at a FitzHugh-Nagumo Solution


X = [x2 x1 x1.^3 x2.*x1.^2];

parfhn = (X'*X)\X'*y;
axes('fontsize',14)
plot3(x1,x2,y-X*parfhn);
xlabel('\fontsize{20} v')
ylabel('\fontsize{20} dv')
zlabel('\fontsize{20} d2v')


axes('fontsize',14)
plot3(x1,x2,y);
hold on
plot3(x1,x2,X*parfhn,'r')
xlabel('\fontsize{20} v')
ylabel('\fontsize{20} dv')
zlabel('\fontsize{20} d2v')



% Now fill this in

X = [o1 x1 x2 x1.*x2 x1.^2 x1.^2.*x2 x1.^3];

parfhn2 = (X'*X)\X'*y;

plot3(x1,x2,y-X*parfhn2);
xlabel('x1')
ylabel('dx1')


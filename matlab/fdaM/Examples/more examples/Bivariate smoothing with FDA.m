%   Bivariate smoothing with FDA

%   Unfortunately the FDA package does not include high-level smoothing
%   functions for bivariate data. 

%   Fortunately, it's not too hard to make it work, anyway. 

%   First, we will load the california income data

load 'calsmall.mat'

whos

%%   These data give income (in $100,000) in terms of latitude and 
%   longitude in California. The first column represents income

hist(calsmall(:,1))

%%   the second are latitude and longitude

plot(calsmall(:,3),calsmall(:,2),'.')

%%   As one might expect, there is a considerable rise in income in cities

plot3(calsmall(:,3),calsmall(:,2),calsmall(:,1),'.')

%%   which makes it tough to be a grad student there. 

%   We will want to smooth these with b-splines for which we need two bases. 

%   First the ranges

max(calsmall)
min(calsmall)

lonrange = [-124.35 -114.31];
latrange = [32.54 41.95];

%%   I am NOT going to use 20640 basis functions in each dimension!

latnknots = 21;
lonnknots = 21;

latnorder = 4;
lonnorder = 4;

latnbasis = latnknots+latnorder-2;
lonnbasis = lonnknots+lonnorder-2;

lonknots = linspace(lonrange(1),lonrange(2),lonnknots);
latknots = linspace(latrange(1),latrange(2),latnknots);

latbasis = create_bspline_basis(latrange,latnbasis,latnorder,latknots);
lonbasis = create_bspline_basis(lonrange,lonnbasis,lonnorder,lonknots);

%%   It is easiest if I use an eliptical penalty, so we'll define

latLfd = int2Lfd(2);
lonLfd = int2Lfd(2);

%%   Now I will calculate the design and penalty matrices associated with each
%   basis

latDmat = eval_basis(calsmall(:,2),latbasis);
lonDmat = eval_basis(calsmall(:,3),lonbasis);

latPmat2 = eval_penalty(latbasis,latLfd);
lonPmat2 = eval_penalty(lonbasis,lonLfd);

latPmat0 = eval_penalty(latbasis,int2Lfd(0));
lonPmat0 = eval_penalty(lonbasis,int2Lfd(0));

%%   I put these together with kronecker products:

Dmat = kron(lonDmat,ones(1,latnbasis)).*kron(ones(1,lonnbasis),latDmat);

D2mat = Dmat'*Dmat;

yvec = Dmat'*calsmall(:,1);

%%   Do a straight smooth of the data

bhat = inv(D2mat)*yvec;

B = reshape(bhat,latnbasis,lonnbasis);

cal_bifd = bifd(B,latbasis,lonbasis);

latfine = linspace(latrange(1),latrange(2),101);
lonfine = linspace(lonrange(1),lonrange(2),101);

bifdvals = eval_bifd(latfine,lonfine,cal_bifd);

surf(lonfine,latfine,bifdvals)
hold on
plot3(calsmall(:,3),calsmall(:,2),calsmall(:,1),'.')

%%   Now define penalty matrices

latPmat = kron(lonPmat0,latPmat2);
lonPmat = kron(lonPmat2,latPmat0);


lambda1 = 0.1;
lambda2 = 0.1;

bhat = inv(D2mat+lambda1*latPmat+lambda2*lonPmat)*yvec;

B = reshape(bhat,latnbasis,lonnbasis);

cal_bifd = bifd(B,latbasis,lonbasis);

latfine = linspace(latrange(1),latrange(2),101);
lonfine = linspace(lonrange(1),lonrange(2),101);

bifdvals = eval_bifd(latfine,lonfine,cal_bifd);

surf(lonfine,latfine,bifdvals)
hold on
plot3(calsmall(:,3),calsmall(:,2),calsmall(:,1),'.')


%%   Perhaps I want to put in some boundaries; here I will just leave out
%   bases for which there is no support. 

leavout = sum(Dmat) == 0;
leavein = 1:size(Dmat,2);
leavein(leavout==1) = [];

D2matl = D2mat(leavein,leavein);
latPmatl = latPmat(leavein,leavein);
lonPmatl = lonPmat(leavein,leavein);

yvecl = yvec(leavein);

bhat = inv(D2matl+lambda1*latPmatl+lambda2*lonPmatl)*yvecl;
bhat2 = zeros(latnbasis*lonnbasis,1);
bhat2(leavein) = bhat;

%%   Now I can go ahead again

B = reshape(bhat2,latnbasis,lonnbasis);

cal_bifd = bifd(B,latbasis,lonbasis);

latfine = linspace(latrange(1),latrange(2),101);
lonfine = linspace(lonrange(1),lonrange(2),101);

bifdvals = eval_bifd(latfine,lonfine,cal_bifd);

surf(lonfine,latfine,bifdvals)
hold on
plot3(calsmall(:,3),calsmall(:,2),calsmall(:,1),'.')

%%   And lets look at some residuals


pred = Dmat*bhat2;
plot3(calsmall(:,3),calsmall(:,2),pred,'.')

hist(calsmall(:,1)-pred)

plot3(calsmall(:,3),calsmall(:,2),calsmall(:,1)-pred,'.')



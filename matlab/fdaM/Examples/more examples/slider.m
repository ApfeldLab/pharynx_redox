load 'slider.mat'
%%

plot(slider)

%%   Now lets define a basis

norder = 5;
knots = 0:4:400;

nbasis = length(knots)+norder - 2;

bbasis = create_bspline_basis([0 400],nbasis,norder,knots);

bfdPar = fdPar(bbasis,int2Lfd(3));

[sfd,df,gcv] = smooth_basis(0:400,slider,bfdPar);


for i = 1:10
    %  i
    bfdPar = fdPar(bbasis,int2Lfd(3),10^(6-i));
    % originql: [sfd,df,gcv] = smooth_basis([0:400],slider,bfdPar,ones(size(slider,1),1),1.4);
    [sfd,df,gcv] = smooth_basis([0:400],slider,bfdPar);
    G(i) = sum(gcv);
    disp(sprintf('i: %.0f, G: %.2f',i,G(i)));
end

bfdPar = fdPar(bbasis,int2Lfd(3),10^3);
[sfd,df,gcv,coef,sse,penmat,y2cmap] = smooth_basis(0:400,slider,bfdPar);


%%   Alright, now we want a linear model

xcell = {pchange,ones(size(pchange))};

betacell = cell(2,1);
betacell(:) = {fdPar(bbasis,int2Lfd(3))};

%[betaestcell,yhatfdobj,betastderrcell,bvariance,c2bmap] = fRegress(sfd,xcell,betacell,y2cmap);
fRegressStruct = fRegress(sfd,xcell,betacell,ones(16,1),y2cmap);

betaestcell = fRegressStruct.betahat;
yhatfdobj   = fRegressStruct.yhat;



%%
figure
plot(getfd(betaestcell{1}))
figure
plot(getfd(betaestcell{2}))
%%

plotfit_fd(moves,0:400,yhatfdobj)
%%

surf(bvariance)
%%

beta1 = eval_fd(0:400,getfd(betaestcell{1}));
beta2 = eval_fd(0:400,getfd(betaestcell{2}));

bstd1 = eval_fd(0:400,betastderrcell{1});
bstd2 = eval_fd(0:400,betastderrcell{2});

plot(beta1)
hold on
plot(beta1+bstd1,'--')
plot(beta1-bstd1,'--')

plot(beta2,'r')
plot(beta2+bstd2,'r--')
plot(beta2-bstd2,'r--')


%%   What happens if I fit pointwise?

svals = eval_fd(0:400,sfd);

x = [pchange ones(size(pchange))];

bhat = inv(x'*x)*(x'*slider');

plot(bhat','m')

hold off

%%   Now lets try a concurrent linear model

dsvals = eval_fd(0:400,sfd,1);

[dsfd,blah,gcv,coef,sse,penmat,y2cmap] = smooth_basis(0:400,dsvals,bfdPar);

plot(dsfd)

%%   Now lets add sfd to the model

xcell = {pchange,sfd};

betacell = cell(2,1);
betacell(:) = {fdPar(bbasis,int2Lfd(3))};

[betaestcell,yhatfdobj,betastderrcell,bvariance,c2bmap] = fRegress(dsfd,xcell,betacell,y2cmap);

plot(getfd(betaestcell{1}))
plot(getfd(betaestcell{2}))

plotfit_fd(dsvals,0:400,yhatfdobj)

surf(bvariance)


beta1 = eval_fd(0:400,getfd(betaestcell{1}));
beta2 = eval_fd(0:400,getfd(betaestcell{2}));

bstd1 = eval_fd(0:400,betastderrcell{1});
bstd2 = eval_fd(0:400,betastderrcell{2});

plot(beta1)
hold on
plot(beta1+bstd1,'--')
plot(beta1-bstd1,'--')

plot(beta2,'r')
plot(beta2+bstd2,'r--')
plot(beta2-bstd2,'r--')


%%   Now lets do some cross validation

scoef = getcoef(sfd);
dscoef = getcoef(dsfd);
lambdas = 10.^(7:12);

for j = 1:length(lambdas)
    
    G = 0;
    for i = 1:size(scoef,2)
        disp([j,i])
        
       tscoef = scoef;
       tscoef(:,i) = [];
       tsfd = fd(tscoef,bbasis);
       
       tdscoef = dscoef;
       tdscoef(:,i) = [];
       tdsfd = fd(tdscoef,bbasis);
       
       tdf = df;
       tdf(i) = [];
       
       xcell = {tdf,tsfd};
       
       betacell = {fdPar(bbasis,int2Lfd(3),lambdas(j)),fdPar(bbasis,int2Lfd(3))};
       betaestcell = fregress(tdsfd,xcell,betacell);
        
       beta1 = eval_fd(0:400,getfd(betaestcell{1}));
       beta2 = eval_fd(0:400,getfd(betaestcell{2}));
       
       yhat = df(i)*beta1 + svals(:,i).*beta2;
       
       G = G + mean( (dsvals(:,i) - yhat).^2 );
       
    end
    
    cv(j) = G;
    
end

cv

%%   It looks like lambda = 10^10 does it

betacell = {fdPar(bbasis,int2Lfd(3),10^11),fdPar(bbasis,int2Lfd(3))};
xcell = {df,sfd};

[betaestcell,yhatfdobj,betastderrcell,bvariance,c2bmap] = fregress(dsfd,xcell,betacell,y2cmap);

plotfit_fd(dsvals,0:400,yhatfdobj)

contour(bvariance)


beta1 = eval_fd(0:400,getfd(betaestcell{1}));
beta2 = eval_fd(0:400,getfd(betaestcell{2}));

bstd1 = eval_fd(0:400,betastderrcell{1});
bstd2 = eval_fd(0:400,betastderrcell{2});

plot(beta1)
hold on
plot(beta1+bstd1,'--')
plot(beta1-bstd1,'--')

plot(beta2,'r')
plot(beta2+bstd2,'r--')
plot(beta2-bstd2,'r--')



%%   Now lets do the principle differential analysis


xfdcell = {sfd};
ufdcell = {fd(pchange',create_constant_basis([0 400]))};

cbasis = create_bspline_basis([0 400],21);

bwtcell = cell(1,2);
bwtcell(:) = {fdPar(cbasis,int2Lfd(2))};

awtcell = cell(1,1);
awtcell(:) = {fdPar(cbasis,int2Lfd(2))};

difeorder = 2;
nfine = 2000;

[bfdcell,afdcell,resfdcell] = ...
    pdacell(xfdcell, bwtcell, awtcell, ufdcell, difeorder,nfine);


%%   In order to iterate, I'll need to pre-define my basis values

bvals = eval_basis(0:400,bbasis);
bmat = bvals'*bvals;

yvals = bvals'*slider;

qvals = 0:0.1:400;

lambda = 10^3;

%%   Now I can set up a new iteration

thisLfd = Lfd(2,bfdcell);
penmat = eval_penalty(bbasis,thisLfd);

Lbvals = eval_basis(qvals,bbasis,thisLfd);
fvals = diag(eval_fd(qvals,getfd(afdcell{1})))*eval_fd(qvals,ufdcell{1});

yfvals = Lbvals'*fvals;

newcvals = inv(bmat + lambda*penmat)*(yvals+lambda*yfvals);

newfd = fd(newcvals,bbasis);

plotfit_fd(moves,0:400,newfd)

%%   Now try the pda again

xfdcell2 = {newfd};

[bfdcell2,afdcell2,resfdcell2] = ...
    pdacell(xfdcell2, bwtcell, awtcell, ufdcell, difeorder,nfine);






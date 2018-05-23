

[t,hmm] = ode45(@runzlefunode,0:0.1:20,[-0.8 0.2],[],wpars);

Tcell{1} = t;
Ycell{1} = hmm(:,1);

range = [0 10];

nbasis = nknots + norder - 2;        
knots = linspace(range(1),range(2),nknots);

quadvals = MakeQuadPoints(knots,nquad);   

basis_obj = MakeBasis(range,nbasis,norder,knots,quadvals,3);
Lfd_obj = fdPar(basis_obj,3,lambda0);

Lfd_cell = {Lfd_obj Lfd_obj};

DEfd1 = smoothfd_cell(Ycell,Tcell,Lfd_cell);

% plot the smooth

figure(2)
axes('fontsize',14)
devals = eval_fdcell(tfine,DEfd1,0);
plot(tfine,devals{1},'r','LineWidth',2);
hold on;
plot(Tcell{1},Ycell{1},'b.');
hold off;



basis_cell = getcellbasis(DEfd1);


startcoefs2 = lsqnonlin(@SplineCoefErr_DEfit,getcoef(DEfd1{2}),[],[],...
    lsopts_out,DEfd1,2,fn,startpars,[]);

DEfd2 = update_fdcell(startcoefs2,2,DEfd1);


newcoefs = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_out,basis_cell,Ycell,Tcell,wts,lambda,fn,[],wpars);

tDEfd = Make_fdcell(newcoefs,basis_cell);

% plot results along with exact solution

figure(3)
devals = eval_fdcell(tfine,tDEfd,0);
axes('fontsize',14)
for i = 1:length(Ycell)
    subplot(length(Ycell),1,i);
    plot(Tcell{i},Ycell{i},'b.');
    hold on
    plot(tfine,devals{i},'r');
    hold off
end

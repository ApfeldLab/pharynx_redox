d2Gdc2 = make_d2gdc2(DEfd,fn,Tcell,wts,lambda,startpars,alg,Afn);

d2Gdcdp = make_d2gdcdp(DEfd,fn,lambda,startpars,alg,Afn,Ycell,Tcell,wts);
    
dcdp = -d2Gdc2\d2Gdcdp;

d2Hdcdp = make_d2SSEdcdp(DEfd,pars,Afn,Ycell,Tcell,wts);

d2Hdc2 = make_d2SSEdc2(DEfd,Tcell,wts,startpars,Afn);

A = Afn.fn(startpars);

[f,J] = djdc_pts(DEfd,Ycell,Tcell,wts,A);
c = J'*f;
E = sum(f.^2);

eps = 1e-6;

tc= 0*c;
tJ = 0*J;
tG = 0*d2Hdc2;

for i = 1:length(coefs)
    tcoef = coefs;
    tcoef(i) = coefs(i) + eps;
    tDEfd = Make_fdcell(tcoef,basis_cell);
    [tf,ttJ] = djdc_pts(tDEfd,Ycell,Tcell,wts,A);
    tc(i) = (sum(tf.^2) - E)/eps;
    tJ(:,i) = (tf-f)/eps;
    tG(:,i) = (J'*tf-c)/eps;
end



devals = [];
for i = 1:length(Ycell)
    n = length(wts{i});
    devals = [devals; -spdiags(sqrt(wts{i}),0,n,n)*...
        cell2mat(eval_fdcell(Tcell{i},DEfd))*squeeze(Afn.dA(i,:,:))];
end


tdevals = 0*devals;

for i = 1:length(startpars)
    tpars = startpars;
    tpars(i) = startpars(i)+eps;
    
    tA = Afn.fn(tpars);
    
    disp( (tA-A)/eps )
    
    tf =  djdc_pts(DEfd,Ycell,Tcell,wts,tA);

    tdevals(:,i) = (tf-f)/eps;    
end




d2H = make_d2SSEdcdp(DEfd,startpars,Afn,Ycell,Tcell,wts);

d2G = 0*d2H;

A = Afn.fn(startpars);
[f,Zmat] = djdc_pts(DEfd,Ycell,Tcell,wts,A);
dgdc = -Zmat'*f;

for i = 1:length(startpars)
    tpars = startpars;
    tpars(i) = startpars(i) + eps;    
    tA = Afn.fn(tpars);
    [terr,tzmat] = djdc_pts(DEfd,Ycell,Tcell,wts,tA);
    
    tdgdc = -tzmat'*terr;
    d2G(:,i) = (tdgdc-dgdc)/eps;
end


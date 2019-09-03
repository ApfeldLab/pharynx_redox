prof_data = squeeze(prof(:,1,1,:))';

norder = 6;
nbasis = length(prof_data) + norder - 2;
heightbasis = create_bspline_basis([1, 200], nbasis, norder);


loglam   = -10:.25:-4;
gcvsave = zeros(length(loglam),1);
dfsave  = gcvsave;

for i=1:length(loglam)
    lambdai   = 10^loglam(i);
    hgtfdPari = fdPar(heightbasis, 4, lambdai);
    [hgtfdi, dfi, gcvi] = smooth_basis(1:200, prof_data, hgtfdPari);
    gcvsave(i) = sum(gcvi);
    dfsave(i)  = dfi;
end

plot(loglam, gcvsave);
scatter(1:200, prof_data(:,1));

%%
% scatter(linspace(1,100,200), prof_data(:,1)); 

i470 = squeeze(prof(:,1,1,:))';
i410 = squeeze(prof(:,2,1,:))';

z410 = (i410 - mean(i410,1)) ./ std(i410,0,1);
z470 = (i470 - mean(i470,1)) ./ std(i470,0,1);

% wormFd = makeWormFd_SJ(prof_data, 'lambda', 10^-10, 'n_breaks', 50);
f410 = makeWormFd_SJ(z410, 'lambda', 10^-.5, 'n_order', 5 ,'n_breaks', 100);
f470 = makeWormFd_SJ(z470, 'lambda', 10^-.5, 'n_order', 5 ,'n_breaks', 100);

figure;
hold on;
scatter(linspace(1,100,200), z410(:,1));
scatter(linspace(1,100,200), z470(:,1));
plot(f410(1));
plot(f470(1));
hold off;

figure;
hold on;
plot(f410(1), 1);
plot(f470(1), 1);
hold off;
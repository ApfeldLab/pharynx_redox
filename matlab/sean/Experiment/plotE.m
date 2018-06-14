exp = e;
% Plot all E
figure;
subplot(2,2,1);
plot(exp.reg.E);

% Plot mean +/- 1.96*std. dev. of E

subplot(2,2,2);
x = linspace(1,100,size(exp.reg.E,1));
% Plot mean for each strain
set(0,'DefaultAxesColorOrder',brewermap(numel(exp.strains),'RdBu'))
for i=1:numel(exp.strains)
    strain = exp.strains{i};
    logIdx = ismember(exp.metadata.Strain, strain);
    strainSubset = exp.reg.E(:,logIdx);

    means = mean(strainSubset.');
    stds = std(strainSubset.');

    boundedline(x, means, 1.96*stds, 'alpha');
    hold all;
end
legend(exp.strains);
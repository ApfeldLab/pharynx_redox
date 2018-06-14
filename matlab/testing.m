% Plot and color according to mean medial_axis value
n_obs = numel(e.reg.regions.E.medial_axis);

map = brewermap(n_obs,'RdYlBu');
figure;
[ma_sorted, ma_order] = sort(e.reg.regions.i470.pm5, 'descend');

sortedE = e.reg.i470(:,ma_order);

axes('ColorOrder',map,'NextPlot','replacechildren') 
for i=1:n_obs
    plot(sortedE(:,i)); 
    hold all
end
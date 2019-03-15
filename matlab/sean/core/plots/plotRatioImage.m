cmap_ = cbrewer('div', 'RdBu', 5, 'PCHIP');
figure;


colorbar;
for i=1:nAnimal1
ratioImage = imR(:,:,i);
flImage = crop_410_1(:,:,i);
R_adj = ja_adjust_brightness(ratioImage, flImage, 4000, cmap_, 1-pd.sigma*3, 1+pd.sigma*3);
image(gca, R_adj);
title(int2str(i));
colormap(map);
caxis([1-pd.sigma*3 1+pd.sigma*3]);
colorbar;
% waitforbuttonpress();
export_fig(sprintf('/Users/sean/Desktop/ratios/%d.pdf', i));
end

%%
map = cbrewer('div', 'RdBu',5,'PCHIP');
x=linspace(1-pd.sigma*2,1+pd.sigma*2,10);
[X,Y] = meshgrid(x);
Z = Y;
figure;
surf(X,Y,Z);
caxis([1-pd.sigma*3 1+pd.sigma*3]);
colormap(map);
colorbar;
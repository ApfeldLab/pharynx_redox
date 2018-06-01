%% Load Data
load handel; % Hallelujah (bell)

ROOT_DIR = '/Users/sean/code/wormAnalysis/data/vab1_2do_05_25_mvmt/';
STRAIN_NAME = 'vab1_2do_05_25_mvmt';

m_410 = loadIntensityData_SJ(strcat(ROOT_DIR, 'PA-vab1_2do_05_25_mvmt-410-measurements.csv'), 0);
m_470 = loadIntensityData_SJ(strcat(ROOT_DIR, 'PA-vab1_2do_05_25_mvmt-470-measurements.csv'), 0);

im_410 = tiffread2(strcat(ROOT_DIR, 'PA-vab1_2do_05_25_mvmt-410.tif'));
im_470 = tiffread2(strcat(ROOT_DIR, 'PA-vab1_2do_05_25_mvmt-470.tif'));

n_worms = size(m_410,2);

coords_410 = flipud(csvread(strcat(ROOT_DIR, 'PA-vab1_2do_05_25_mvmt-410-polyCoords.csv'), 1, 1));
coords_410_x = coords_410(:, 1:2:end);
coords_410_y = coords_410(:, 2:2:end);

coords_470 = flipud(csvread(strcat(ROOT_DIR, 'PA-vab1_2do_05_25_mvmt-470-polyCoords.csv'), 1, 1));
coords_470_x = coords_470(:, 1:2:end);
coords_470_y = coords_470(:, 2:2:end);

%% Register
wormFds = cell(n_worms,1);
regWormFds = cell(n_worms,1);
warpFds = cell(n_worms,1);

WARP_NBASIS = 6;
WARP_ORDER = 4;
WARP_LAMBDA = 5000;


parfor i=1:n_worms
    wormFds(i) = {makeWormFd_SJ(horzcat(m_410(:,i),m_470(:,i)), 'lambda', 10^0.0891)};
    warpBasis = create_bspline_basis([1 100], WARP_NBASIS, WARP_ORDER);
    [cregFd, cwarpFd] = register_fd(wormFds{i}(1), wormFds{i}, fdPar(warpBasis, int2Lfd(2), WARP_LAMBDA));
    regWormFds(i) = {cregFd};
    warpFds(i) = {cwarpFd};
end

%% Plot
xs = linspace(1,100,1000);
rows = 3;
cols = 2;

f = figure('units','normalized','outerposition', [0 0 1 1]);
set(f, 'Visible', 'off');

SAVE_FIGS = 1;

for i=20:23
    subplot(rows,cols,1);
    
    plot(eval_fd(xs, wormFds{i}(1)), 'Color', 'b'); hold on;
    plot(eval_fd(xs, wormFds{i}(2)), 'Color', 'r'); hold on;
    plot(eval_fd(xs, regWormFds{i}(2)), 'r-.'); hold off;
    title(strcat(num2str(i)));
    legend('410', '470' , 'Reg470', 'Location', 'northeast');
    
    subplot(rows,cols,2);
    plot(eval_fd(xs, warpFds{i}));
    title(strcat('B(', num2str(WARP_NBASIS),'); ',...
        'O(', num2str(WARP_ORDER), '); ',...
        'L(', num2str(WARP_LAMBDA),')'));
    
%     subplot(rows,cols,2);
%     plot(eval_fd(xs, regWormFds{i}));
%     title(strcat('Registered (', num2str(i), ')'));
    
    subplot(rows,cols,3);
    imagesc(im_410(i).data); hold on;
    plot(coords_410_x(:,i), coords_410_y(:,i), 'LineWidth', 3, 'Color', 'b');
    title('410');
    
    subplot(rows,cols,4);
    imagesc(im_470(i).data); hold on;
    plot(coords_470_x(:,i), coords_470_y(:,i), 'LineWidth', 3, 'Color', 'r');
    title('470');
    
    subplot(rows,cols,5);
    imagesc(double(im_410(i).data)./double(im_470(i).data)); hold on;
    plot(coords_410_x(:,i), coords_410_y(:,i), 'LineWidth', 3, 'Color', 'b'); hold on;
    plot(coords_470_x(:,i), coords_470_y(:,i), 'LineWidth', 3, 'Color', 'r');
    legend('410', '470', 'Location', 'northwest');
    title('410/470');
    colormap(rgb);
    colorbar
    set(gca,'CLim',[.5 1.5]);
    
    subplot(rows,cols,6);
    plot(eval_fd(xs, wormFds{i}(1)), 'Color', 'b'); hold on;
    plot(eval_fd(xs, regWormFds{i}(2)), 'Color', 'r'); hold off;
    title(strcat(num2str(i)));
    legend('410', 'Reg470', 'Location', 'northeast');
%     suplabel('Single title on top', 't');
    if SAVE_FIGS == 1 
        export_fig(sprintf(strcat(ROOT_DIR, 'figs/%d.pdf'), i));
    end
    clf('reset');
end

% %%
% figure;
% 
% R = zeros(n_worms,1000);
% for i=1:n_worms
%     rFd = regWormFds{i}(1) ./ regWormFds{i}(2);
%     R(i,:) = eval_fd(xs, rFd);
% end
% 
% plot(xs, mean(ja_E(ja_oxd(R))));
%% Load Data
load handel; % Hallelujah (bell)

m_410 = loadIntensityData_SJ('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_410_measure.csv', 0);
m_470 = loadIntensityData_SJ('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_470_measure.csv', 0);

im_410 = tiffread2('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_410_PA.tif');
im_470 = tiffread2('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_470_PA.tif');

n_worms = size(m_410,2);

coords_410 = flipud(csvread('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_410_polycoords.csv', 1, 1));
coords_410_x = coords_410(:, 1:2:end);
coords_410_y = coords_410(:, 2:2:end);

coords_470 = flipud(csvread('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_470_polycoords.csv', 1, 1));
coords_470_x = coords_470(:, 1:2:end);
coords_470_y = coords_470(:, 2:2:end);

%% Register
[fdObjs, data] = ChannelRegister(ssquare(m_410), ssquare(m_470));

%%
thresh = 1500;
data.m410 = ssquare(clip_sj(data.m410, thresh));
data.m470 = ssquare(clip_sj(data.m470, thresh));

plot(data.m470);
plotRegionBoundaries(Constants.regions);

%% Plot
rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;
b = repmat(linspace(0,1,200),20,1);

xs = linspace(1,100,1000);
rows = 3;
cols = 2;

plots = 0;
if plots == 1
    f = figure('units','normalized','outerposition', [0 0 1 1]);
    set(f, 'Visible', 'on');
    for i=1:n_worms
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
    colorbar;
    set(gca,'CLim',[.5 1.5]);
    
    subplot(rows,cols,6);
    plot(eval_fd(xs, wormFds{i}(1)), 'Color', 'b'); hold on;
    plot(eval_fd(xs, regWormFds{i}(2)), 'Color', 'r'); hold off;
    title(strcat(num2str(i)));
    legend('410', 'Reg470', 'Location', 'northeast');
    pause;
%     suplabel('Single title on top', 't');
%     export_fig(sprintf('sean/figs/channel_registration/SAY98_eat5_2do_05_24_18/%d.pdf', i));
    end
end

sound(y,Fs); % Ring the bell!
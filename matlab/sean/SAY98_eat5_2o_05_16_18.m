%% Load Data
m_410 = loadIntensityData_SJ('../data/SAY98_eat5_2do_05_16/SAY98_eat5_2do_05_16_410.csv', 0);
m_470 = loadIntensityData_SJ('../data/SAY98_eat5_2do_05_16/SAY98_eat5_2do_05_16_470.csv', 0);

im_410 = tiffread2('../data/SAY98_eat5_2do_05_16/PA-SAY98_eat5_2d0_2018_05_16_410.tif');
im_470 = tiffread2('../data/SAY98_eat5_2do_05_16/PA-SAY98_eat5_2d0_2018_05_16_470.tif');

n_worms = size(m_410,2);

%% Register
wormFds = cell(n_worms,1);
regWormFds = cell(n_worms,1);

parfor i=1:40
    wormFds(i) = {makeWormFd_SJ(horzcat(m_410(:,i),m_470(:,i)), 'lambda', 0.0891)};
%     [cregFd, cwarpFd] = register_fd(wormFds{i}(1), wormFds{i}, fdPar(getbasis(wormFds{i}), int2Lfd(2), 1000));
    regWormFds(i) = {wormFds(i)};
end

%% Colors
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

%% Plot
figure;
xs = linspace(1,100,1000);
rows = 3;
cols = 2;

for i=1:49
    subplot(rows,cols,1);
    plot(eval_fd(xs, wormFds{i}));
    title(strcat('Unregistered (', num2str(i), ')'));
    legend('410', '470');
    
%     subplot(rows,cols,2);
%     plot(eval_fd(xs, regWormFds{i}));
%     title(strcat('Registered (', num2str(i), ')'));
    
    subplot(rows,cols,3);
    imagesc(fliplr(im_410(i).data));
    title('410');
    
    subplot(rows,cols,4);
    imagesc(fliplr(im_470(i).data));
    title('470');
    
    subplot(rows,cols,5);
    imagesc(fliplr(double(im_410(i).data)./double(im_470(i).data)));
    title('410/470');
    colormap(rgb);
    colorbar
    set(gca,'CLim',[.5 1.5])
    pause;
end
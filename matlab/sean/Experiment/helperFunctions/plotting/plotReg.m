figures = [];

% Generate Figures
for i=1:1
%     figures(i) = figure('visible','off');
    figures(i) = figure;
    % Mean i410
    ax = subplot(3,3,1);
    e.plotMeanByStrain('reg', 'i410', ax);
    title('Mean i410 +/- 1.96*std');

    % Mean i470
    ax = subplot(3,3,2);
    e.plotMeanByStrain('reg', 'i470', ax);
    title('Mean i470 +/- 1.96*std');

    % Mean E
    ax = subplot(3,3,3);
    e.plotMeanByStrain('reg', 'E', ax);
    title('Mean E +/- 1.96*std');

    % im410/im470
    % subplot(3,3,4);
    % imagesc(im470(i).data ./ im410(i).data);

    % Raw vs. Registered 470
    xs = linspace(1,100,1000);
    subplot(3,3,5);
    plot(1:100, e.raw.sq410(:,i), 1:100, e.raw.sq470(:,i), xs, e.reg.i470(:,i));
    legend('raw 410', 'raw470', 'reg470');
    title(sprintf('Raw vs. Registered 470 (%d)', i));

    % Registered 410 & mean of Strain
    ax = subplot(3,3,7);
    strain = e.metadata.Strain{i};
    e.plotMeanByStrain('reg', 'i410', ax); hold on;
    plot(xs, e.reg.i410(:,i), '--','DisplayName',strcat('Reg. Worm', num2str(i))); hold off;
    title(sprintf('Mean 410nm +/- 1.96*std of %s and worm #%d', strain, i));

    % Registered 470 & mean of Strain
    ax = subplot(3,3,8);
    strain = e.metadata.Strain{i};
    e.plotMeanByStrain('reg', 'i470', ax); hold on;
    plot(xs, e.reg.i470(:,i), '--','DisplayName',strcat('Reg. Worm', num2str(i))); hold off;
    title(sprintf('Mean 470nm +/- 1.96*std of %s and worm #%d', strain, i));

    % Registered E & mean of Strain
    ax = subplot(3,3,9);
    strain = e.metadata.Strain{i};
    e.plotMeanByStrain('reg', 'E', ax); hold on;
    plot(xs, e.reg.E(:,i), '--','DisplayName',strcat('Reg. Worm', num2str(i))); hold off;
    title(sprintf('Mean E +/- 1.96*std of %s and worm #%d', strain, i));
    
    set(gcf, 'Position', get(0, 'Screensize'));
    fname = fullfile(experiment_directory, 'figs', 'registration', strcat(num2str(i), '.pdf'));
    export_fig fname;
    close(gcf);
end


figSize = [21, 29];      % [width, height]
figUnits = 'Centimeters';

%% Save Figures
for i=1:numel(figures)
    fig = figures(i);
     
    % Resize the figure
%     set(fig, 'Units', figUnits);
%     pos = get(fig, 'Position');
%     pos = [pos(1), pos(4)+figSize(2), pos(3)+figSize(1), pos(4)];
%     set(fig, 'Position', pos);
    set(fig, 'Position', [1 41 1920 964]);
    
    fig_name = fullfile(e.directory, 'figs', sprintf('%d.pdf', i));
    saveas(fig,fig_name);
end
function plotChannelRegistration(wormFd ,regFd, warpFd, im410, im470, coords410, coords470)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    figure;
    rows = 3;
    cols = 2;
    
    xs = linspace(1,100,1000);
    
    subplot(rows,cols,1);
    plot(eval_fd(xs, wormFd(1)), 'Color', 'b'); hold on;
    plot(eval_fd(xs, wormFd(2)), 'Color', 'r'); hold on;
    plot(eval_fd(xs, regWormFd(2)), 'r-.'); hold off;
    legend('410', '470' , 'Reg470', 'Location', 'northeast');
    
    subplot(rows,cols,2);
    plot(eval_fd(xs, warpFd));
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
    pause;
%     suplabel('Single title on top', 't');
%     export_fig(sprintf('sean/figs/channel_registration/SAY98_eat5_2do_05_24_18/%d.pdf', i));
end

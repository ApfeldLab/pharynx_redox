plot(xs410, ys410); hold on;
plot(warped_xs470, warped_ys470); hold on;

midpoint = (max(ys410) + min(ys410)) / 2;

downsample_rate = 100;
x = downsample(xs410, downsample_rate);
y = downsample((ones(size(xs410)) * midpoint), downsample_rate);
u = downsample(warped_xs470 - xs410, downsample_rate);
v = downsample(warped_ys470.' - ys410.', downsample_rate);

radians = atan(y ./ x);

quiver(x, y, u, v, 'LineWidth', 2); hold off;
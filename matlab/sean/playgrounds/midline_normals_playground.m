% rebasedVecs = cell(length(norms1), 1);
% cob = cell(length(norms1), points);
% for i = 1:length(norms1)
%     for j = 1:points
%         cob{i,j} = [norms1{i}(j,:).' tangs1{i}(j,:).'];
%         rebasedVecs{i}(j,:) = (cob{i,j}\matchingVecs{i}(j,:).').';
%     end
% end
%% 


[dx, dy] = getDyDxMidlines(mids1, mids2, scaled_bounds1, scaled_bounds2, fdObjs, 100);

%%
[tangs1, norms1] = getNormsAndTangs(mids1, 100);

animal = 21;
xdata = (1:100).';
ydata = feval(mids1{animal}, xdata);
tangs = tangs1{animal};
norms = norms1{animal};


subplot(1,3,1);
hold on;
% plot(xdata, ydata);
quiver(xdata, ydata, norms(:,1), norms(:,2), 'AutoScaleFactor', 0.1);
quiver(xdata, ydata, tangs(:,1), tangs(:,2), 'AutoScaleFactor', 0.1);
quiver(xdata, ydata, matchingVecs{animal}(:,1), matchingVecs{animal}(:,2), 'Color', 'black');
xlim([1 100]);
ylim([1 100]);
pbaspect([1 1 1]);
subplot(1,3,2);
quiver(xdata, zeros(size(xdata)), dx(:,animal), dy(:,animal), 'ShowArrowHead', 'off');
xlim([-4 105]);
pbaspect([1 1 1]);
subplot(1,3,3);
hold on;
plot(dx(:,animal));
plot(dy(:,animal));
hline(0);
legend('dx', 'dy');
pbaspect([1 1 1]);

%%
function [dx, dy] = getDyDxMidlines(midlines1, midlines2, scaledBounds1, scaledBounds2, fdObj, nVecs)
    nAnimals = length(midlines1);
    dx = zeros(nVecs, nAnimals);
    dy = zeros(nVecs, nAnimals);

    matchingVecs = getMidlineMatchingVectors(midlines1, midlines2, scaledBounds1, scaledBounds2, fdObj, nVecs);
    [tangs, norms] = getNormsAndTangs(midlines1, nVecs);
    cob = cell(length(norms), nVecs);
    for i = 1:nAnimals
        for j = 1:nVecs
            cob{i,j} = [norms{i}(j,:).' tangs{i}(j,:).'];
            dxdy = (cob{i,j}\matchingVecs{i}(j,:).').';
            dx(j,i) = dxdy(2);
            dy(j,i) = dxdy(1);
        end
    end
end
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
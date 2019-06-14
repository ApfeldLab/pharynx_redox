function [tangs, norms] = getNormsAndTangs(midlines, nVecs)
    nMids = length(midlines);
    tangs = cell(nMids,1);
    norms = cell(nMids,1);
    
    for i=1:nMids
        m = midlines{i};
        d = differentiate(m, linspace(1, 100, nVecs));
        
        theta = atan(d);
        tangs{i} = [cos(theta) sin(theta)];
        norms{i} = [cos(theta + pi/2) sin(theta + pi/2)];
    end
end
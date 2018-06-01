function coords = loadCoordinates(filename)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    allcoords = dlmread(filename, '', 1, 1);
    coords = struct;
    coords.x = ssquare(allcoords(:,1:2:end));
    coords.y = ssquare(allcoords(:,2:2:end));
end
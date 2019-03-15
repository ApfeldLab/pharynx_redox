function coords = loadCoordinates(filename)
%LOADCOORDINATES Load the coordinates of the midlines from FIJI.
%   Returns a struct where the x attribute contains all x-coordinates, and
%   the y attribute contains all y-coorinates. In each, every COLUMN is a
%   different animal.

    allcoords = readtable(filename, 'ReadVariableNames', 0, 'HeaderLines', 1);
    
    % First column is an index... this comes from FIJI. Get rid of it here.
    allcoords = allcoords(:, 2:end);
    
    coords = struct;
    coords.x = table2array(allcoords(:,1:2:end));
    coords.y = table2array(allcoords(:,2:2:end));
end
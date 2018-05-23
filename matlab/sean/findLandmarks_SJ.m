function landmarks = findLandmarks_SJ(fdObj)
%findLandmarks Return x-position of landmarks for each column of intensity 
%data for use in function registration
%   Landmarks:
%       - Maximum in first half of vector
%       - Maximum in second half of vector

basis_range = getbasisrange(getbasis(fdObj));
xvec = linspace(basis_range(1), basis_range(2), 1000);
yvec = eval_fd(xvec, fdObj);

[~, max1] = max(yvec(1:end/2,:));
[~, max2] = max(yvec(end/2:end,:));
max2 = max2 + size(yvec, 1) / 2;

landmarks = vertcat(xvec(max1), xvec(max2)).';
end
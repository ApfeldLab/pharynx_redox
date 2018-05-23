function coefs = getcellcoefs(fd_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getcellcoefs
%
% Returns a vector obainted by concatenating the coefficients taken from
% the functional data objects contained in the cell array fd_cell. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fd_cell = fd_cell';

coefs = getcoef(fd_cell{1});

for i = 2:numel(fd_cell) 
    coefs = [coefs; getcoef(fd_cell{i})];
end

end
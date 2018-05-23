function M = cell2mat(H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function cell2mat
%
% Takes a cell array H of matrices and concatenates them into one matrix
%
% Assumes that the dimensions of the matrices allow this in the sense that
% matrices in the same row of H have the same number of rows and matrices
% in the same column of H have the same number of columns. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = [];

M2 = cell(size(H,2),1);

for i = 1:size(H,1)
    M2{i} = [];
    for j = 1:size(H,2)
        M2{i} = [M2{i},H{i,j}];
    end
    M = [M; M2{i}];
end

end
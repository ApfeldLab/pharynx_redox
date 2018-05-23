function r = forcingdfdp(t,fd_cell,p,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingdfdp
%
% forcing function estimation
%
% derivative with respect to forcing coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = eval_basis_cell(t,more.basisp,0);

r = cell(size(fd_cell,2),length(p));
r(:) = {0};

ind = 1;

for i = 1:length(more.which)
    for j = 1:size(fs{i},2)
        r{more.which(i),ind} = fs{i}(:,j);
        ind = ind+1;
    end
end

end
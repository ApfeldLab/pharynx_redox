function fd_cell = smoothfd_cell(Ycell,Tcell,Lfd_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function smoothfd_cell
%
% Creates a cell array of functional data objects of smooths to a matrix
% of values
%
% INPUTS:
%
% Ycell   -  cell-array or matrix of observations, columns (or cells) 
%           represent components.
%
% Tcell  -  observation times for the obs in Ycell. May be a vector 
%           or cell array if times are different for different components
%
% Lfd_cell - a cell array of Lfd objects; one for each column in Ycell.
%
% OUTPUT
% 
% fd_cell  -  a cell array of functional data objects one for each column
%             of Ycell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Lfd_cell) == 1
    Lfd_cell = repmat(Lfd_cell,size(Ycell,1),size(Ycell,2));
end

fd_cell = cell(size(Ycell));

t2 = Tcell;

for i = 1:numel(Ycell)
    
    if iscell(Tcell)
        t2 = Tcell{i};
    end
    
    if ~isempty(t2)
        fd_obj = smooth_basis(t2,Ycell{i},Lfd_cell{i});
    else
        fd_obj = fd(zeros(getnbasis(getbasis(getfd(Lfd_cell{i}))),1),...
            getbasis(getfd(Lfd_cell{i})));
    end
    fdnames{1} = 'time';
    fdnames{2} = 'dimension';
    fdnames{3} = 'value';

    fd_obj = putnames(fd_obj,fdnames);

    fd_cell{i} = fd_obj;
end

end
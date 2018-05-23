function DEfd_cell = data2fd_cell(basis_cell,Ycell,Tcell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function DEfd_cell
% 
% Creates a cell array of functional data objects to fit data given in a
% cell array Ycell. 
%
% INPUTS:
%
% basis_cell - a cell array of basis objects, one for each component
%
% Ycell - a matrix of observations, each element corresponds to a component
%         of basis_cell
%
% Tcell - times at which the system is observed. May be a vector or cell
%         array if times are different for different components.
%
% OUTPUTS:
%
% DEfd_cell - a cell array of unpenalized fits to the data, each component
%             is a functional data object corresponding to the entries
%             given in Ycell and Tcell. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t2 = Tcell;

DEfd_cell = cell(size(Ycell));

for i = 1:numel(Ycell) 
    if iscell(Tcell) 
        t2 = Tcell{i};
    end
    
    DEfd = data2fd(Ycell{i},t2,basis_cell{i});

    DE_fdnames{1} = 'time';
    DE_fdnames{2} = 'dimension';
    DE_fdnames{3} = 'value';

    DEfd = putnames(DEfd,DE_fdnames);

    DEfd_cell{i} = DEfd;
end

end
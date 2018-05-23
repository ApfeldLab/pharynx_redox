function more = checkmore(more,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkmore
%
% checks the struct more for use with the genlin package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(more)
    more.sub = [kron((1:n)',ones(n,1)) kron(ones(n,1),(1:n)')];
    more.force = [];
    more.mat = zeros(n);
else
    if ~isfield(more,'sub')
        more.sub = [kron((1:n)',ones(n,1)) kron(ones(n,1),(1:n)')];
    end
    if ~isfield(more,'force') more.force = []; end
    if ~isfield(more,'mat') more.mat = zeros(n); end
    
    if ~isempty(more.force)
        m = length(more.force);
        if ~isfield(more,'force_mat') more.force_mat = zeros(n,m); end
        if ~isfield(more,'force_sub')
            more.force_sub = [kron((1:n)',ones(m,1)) kron(ones(n,1),(1:m)')];
        end
        if size(more.force_sub,2)==1
            more.force_sub = [more.force_sub more.force_sub];
        end
        if ~isfield(more,'force_input') more.force_input = []; end

        %        if ~isfield(more,'dforcedp') more.dforcedp = []; end
        %        if ~isfield(more,'d2forcdp2') more.d2forcedp2 = []; end
    end

end
function S = make_sigma(fd_cell,Tcell,Ycell,ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_sigma
%
% Creates a diagonal matrix giving the estimated variance of each
% observation in Ycell. Different components of the Ycell may be permitted
% to have different variances. 
%
% INPUTS:
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% Ycell       - cell array of observed values of the system.
%
% ind        - an indicator function of common variances can take values:
%                0:  each component of each replicate can be different.
%                1:  components accross replicates have the same variance
%                2:  different replicates have different variances, but
%                    components are the same.
%                3:  all comonents have the same variance. 
%
% OUTPUT: 
%
% s      - a diagonal matrix of the size of the number of observations
%           giving the estimated variance/co-variance matrix of the
%           observational errors assuming independence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = zeros(size(fd_cell));
SSE = zeros(size(fd_cell));

S = cell(size(fd_cell));

devals = eval_fdcell(Tcell,fd_cell,0);

for i = 1:size(Ycell,1)
    for j = 1:size(Ycell,2)
        n(i,j) = length(Tcell{i,j});
        if n(i,j)>1
            SSE(i,j) = sum( (Ycell{i,j}-devals{i,j}).^2 );
            if ind == 0
                S{i,j} = SSE(i,j)/(n(i,j)-1)*speye(n(i,j),n(i,j));
            end
        else
            if n(i,j) == 1
                SSE(i,j) = sum( (Ycell{i,j}-devals{i,j}).^2 );
                if ind == 0
                    S{i,j} = SSE(i,j);
                end
            else
                SSE(i,j) = 0;
                if ind == 0
                    S{i,j} = [];
                end
            end
        end 
    end
end

if ind == 1
    n = sum(n,1);
    SSE = sum(SSE,1);
elseif ind==2
        n = sum(n,2);
        SSE = sum(SSE,2);
else
    n = sum(sum(n));
    SSE = sum(sum(SSE));
end

if ind ~= 0
    for i = 1:size(n,1)
        for j = 1:size(n,2)
            S{i,j} = SSE(i,j)/(n(i,j)-1)*speye(n(i,j),n(i,j));
        end
    end
end

S = mattdiag_cell(reshape(S',numel(S),1),0);

end

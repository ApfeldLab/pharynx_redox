function W = SparseJMfun(Jinfo,Y,flag,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SparseJMfun
%
% A sparse matrix multiplication routine (actually, just plain matrix
% multiplication), for use in lsqnonlin.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = Jinfo;

if flag == 0
    W = J'*(J*Y);
else if flag > 0
        W = J*Y;
    else if flag < 0
            W = J'*Y;
        end
    end
end

end
function r = fhnd2fdx2dp(t,fd_cell,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fhnd2fdy2dp
%
% The third derivative of the FitzHugh Nagumo equations with respect 
% to the input x (twice) and parameters p at time t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = cell([2,2,2,3]);
r(1:2,1:2,1:2,1:3) = {0};
r{1,1,1,3} = -2*eval_fd(t,fd_cell{1});
end
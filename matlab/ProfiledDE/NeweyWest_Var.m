function C = NeweyWest_Var(H,g,maxlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NeweyWest_Var
%
% Creates a Newey-West estimate of variance for profiled estimation
% routines. 
%
% Inputs
%
% H - estimated second derivative of the profile objective function
%
% g - estimated gradient of the profile objective function
%
% maxlag - number of lags to use
%
% Output
% 
% C - estimated covariance up to the scale of the residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 3
        n = size(g,1);
        maxlag = max(5,n^(0.2));
    end
    disp(maxlag)
 %   V = inv(H);
V = H;
    I = 0*H;

    if maxlag > 0
        for i = 1:(size(g,2)-1)
            for j = (i+1):size(g,2)
                I(i,j) = Newey_West(g(:,i),g(:,j),maxlag);
                I(j,i) = I(i,j);
            end
        end
    end

    C = V*(I + g'*g)*V;
end



function out = Newey_West(x,y,maxlag)
%% Basic Newey-West Functionality
    w = 1 - (1:maxlag)/maxlag+1;
    w = w/length(x);
    out = mean(x.*y);
    
    for i = 1:maxlag
       out = out + w(i)*sum(trimr(x,i,0).*trimr(y,0,i))+w(i)*sum(trimr(y,i,0).*trimr(x,0,i));
    end
    
end


function out = trimr(a,n1,n2)
    da = size(a);
    
    if  da(2) == 1 | da(1) == 1, 
        out = a( (n1+1):(length(a)-n2) );
    else 
        out = a( (n1+1):(length(a)-n2),:); 
    end
        
end
function SSE = ApproxForwardCV(pars,DEfd,fn,Y,t,dims,which,extras)

if nargin < 7, which = 1:length(t); end
if nargin < 6, dims = 1:size(Y,2); end

SSE = 0;

if size(which,1) == 1; which = which'; end
if size(which,2) == 1
    which = [which (which+1)];
end

devals = cell2mat( eval_fdcell(t,DEfd) );
ddevals = cell2mat( eval_fdcell(t,DEfd,1) );
fdevals = cell2mat(fn.fn(t,DEfd,pars,extras) ); 

h = t(which(:,2)) - t(which(:,1));

SSE = sum(mean( (Y(which(:,2)) - devals(which(:,2),dims)).^2 + ...
    h.^2.*(ddevals(which(:,1),dims)-fdevals(which(:,1),dims)).^2 + ...
    2*h.*(Y(which(:,2)) - devals(which(:,2),dims)).*(ddevals(which(:,1),dims)-fdevals(which(:,1),dims)) ));


% subplot(10,1,i+1)%one subplot per year
% plot(t,y(:,3))    
% r=[];   
% disp(y)
%     for j=1:52     
%     rn = (y(j,3)-Yvec(j))^2; %squared difference between data and ode simulation for every week of a given year
%     r=[r rn]
%     ss=sum(r)
%     end
% S=[S ss]
% end
% F=sum(S)
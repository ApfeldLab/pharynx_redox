% S=[]
% for i=0:2 % years
% p = newpars; % estimated parameters for function beta.
% m0=cell2mat(devals);
% Yvec=cell2mat(Ycell);
% r0 = (m0(1+52*i,:))'; % initial condition of estimated susceptibles, exposed and infected, at the beginning of each year 
% t0 = 1939+i;
% tf = 1940+i;
% tvec = [t0:1/52:tf];
% options = [];

function SSE = ForwardCV(pars,odefn,devals,Y,t,dims,which,options,extras)

if nargin < 9, extras = []; end
if nargin < 8, options = []; end
if nargin < 7, which = 1:length(t); end
if nargin < 6, dims = 1:size(Y,2); end

SSE = 0;

if size(which,1) == 1; which = which'; end

if iscell(devals)
    if isa_fd(devals{1}) 
        devals = cell2mat( eval_fdcell(t,devals) ); 
    else
        devals = cell2mat(devals);
    end
end

for i = 1:size(which,1)
    if iscell(which)
        tvec = t(which{i});
        Ycomp = Y(which{i},:);
        r0 = devals(which{i}(1),:);
    elseif size(which) > 1
        if size(which,2) == 2
            tvec = t(which(i,1):which(i,2));
            Ycomp = Y(which(i,1):which(i,2),:);
        else
            tvec = t(which(i,:));
            Ycomp = Y(which(i,:),:);
        end
        r0 = devals(which(i,1),:);
    end

    remove = [0; diff(tvec) <= 0];
    Ycomp = Ycomp(~remove,:);
    tvec = tvec(~remove,:);
    
    [tt,ty]=ode23(odefn,tvec,r0',options,pars,extras);
    
    if length(tvec) == 2, 
        ty = ty(size(ty,1),:); 
        Ycomp = Ycomp(2,:);
    end
    
    SSE = SSE + mean(mean( (ty(:,dims)-Ycomp).^2 ));
end    
    
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
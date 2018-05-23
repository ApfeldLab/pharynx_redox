% Tests of derivatives

eps = 1e-6;

% Testing dfdx

DFDX = fn.dfdx(tspan,DEfd,pars,fn.extras);

F = fn.fn(tspan,DEfd,pars,fn.extras);

tDFDX = DFDX;
for i = 1:length(DEfd)
   coefs = getcoef(DEfd{i});
   tDEfd = DEfd;
   tDEfd{i} = fd(coefs+eps,getbasis(DEfd{i}));
   
   tF = fn.fn(tspan,tDEfd,pars,fn.extras);
   for j = 1:length(DEfd)
       tDFDX{j,i} = (tF{j} - F{j})/eps;
   end
    
end


for i = 1:size(DFDX,1)
    for j = 1:size(DFDX,2)

        subplot(2,1,1)
        plot(DFDX{i,j})
        hold on
        plot(tDFDX{i,j},'r')
        hold off
        title(strcat(num2str(i),' ',num2str(j)),'fontsize',20)
        
        subplot(2,1,2)
        plot(tDFDX{i,j} - DFDX{i,j})
        
        pause
        
    end
end



% Now lets try some parameters

DFDP = fn.dfdp(tspan,DEfd,pars,fn.extras);

F = fn.fn(tspan,DEfd,pars,fn.extras);

tDFDP = DFDP;
for i = 1:length(pars)
    tpars = pars;
    tpars(i) = pars(i) + eps;
    
   tF = fn.fn(tspan,DEfd,tpars,fn.extras);
   for j = 1:length(DEfd)
       tDFDP{j,i} = (tF{j} - F{j})/eps;
   end
    
end


for i = 1:size(DFDP,1)
    for j = 1:size(DFDP,2)

        subplot(2,1,1)
        plot(DFDP{i,j})
        hold on
        plot(tDFDP{i,j},'r')
        hold off
        title(strcat(num2str(i),' ',num2str(j)),'fontsize',20)
        
        subplot(2,1,2)
        plot(tDFDP{i,j} - DFDP{i,j})
        
        pause
        
    end
end



% Ok, second derivatives wrt X


D2FDX2 = fn.d2fdx2(tspan,DEfd,pars,fn.extras);

F = fn.dfdx(tspan,DEfd,pars,fn.extras);

tD2FDX2 = D2FDX2;
for i = 1:length(DEfd)
   coefs = getcoef(DEfd{i});
   tDEfd = DEfd;
   tDEfd{i} = fd(coefs+eps,getbasis(DEfd{i}));
   
   tF = fn.dfdx(tspan,tDEfd,pars,fn.extras);
   for j = 1:length(DEfd)
       for k = 1:length(DEfd)
           tD2FDX2{j,k,i} = (tF{j,k} - F{j,k})/eps;
       end
   end 
end


for i = 1:size(D2FDX2,1)
    for j = 1:size(D2FDX2,2)
        for k = 1:size(D2FDX2,3)

            subplot(2,1,1)
            plot(D2FDX2{i,j,k})
            hold on
            plot(tD2FDX2{i,j,k},'r')
            hold off
            title(strcat(num2str(i),num2str(j),num2str(k)),'fontsize',20)

            subplot(2,1,2)
            plot(tD2FDX2{i,j,k} - D2FDX2{i,j,k})

            pause

        end
    end
end


% Ok, second derivatives wrt P


D2FDXDP = fn.d2fdxdp(tspan,DEfd,pars,fn.extras);

F = fn.dfdx(tspan,DEfd,pars,fn.extras);

tD2FDXDP = D2FDXDP;
for i = 1:length(pars)
   tpars = pars;
   tpars(i) = pars(i) + eps;
   
   tF = fn.dfdx(tspan,DEfd,tpars,fn.extras);
   for j = 1:length(DEfd)
       for k = 1:length(DEfd)
           tD2FDXDP{j,k,i} = (tF{j,k} - F{j,k})/eps;
       end
   end 
end


Fp = fn.dfdp(tspan,DEfd,pars,fn.extras);
ttD2FDXDP = D2FDXDP;
for i = 1:length(DEfd)
   coefs = getcoef(DEfd{i});
   tDEfd = DEfd;
   tDEfd{i} = fd(coefs+eps,getbasis(DEfd{i}));
   
   tF = fn.dfdp(tspan,tDEfd,pars,fn.extras); 
    
   for j = 1:length(DEfd)
       for k = 1:length(pars)
           ttD2FDXDP{j,i,k} = (tF{j,k} - Fp{j,k})/eps;
       end
   end 
end




for i = 1:size(D2FDXDP,1)
    for j = 1:size(D2FDXDP,2)
        for k = 1:size(D2FDXDP,3)
            subplot(2,1,1)
            plot(D2FDXDP{i,j,k})
            hold on
            plot(tD2FDXDP{i,j,k},'r')
            hold off
            title(strcat(num2str(i),num2str(j),num2str(k)),'fontsize',20)

            subplot(2,1,2)
            plot(tD2FDXDP{i,j,k} - D2FDXDP{i,j,k})

            pause

        end
    end
end


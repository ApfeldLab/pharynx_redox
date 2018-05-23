
addpath('../lincub')

% I assume we've run FHN_model_building.m

fn.fn       = @lincubfun;       % RHS function

fn.dfdx     = @lincubdfdx;      % Derivative wrt inputs (Jacobian)
fn.dfdp     = @lincubdfdp;      % Derviative wrt parameters

fn.d2fdx2   = @lincubd2fdx2;    % Hessian wrt inputs
fn.d2fdxdp  = @lincubd2fdxdp;   % Hessian wrt inputs and parameters

                          
% parameters

startpars = [newpars; 0];
                             

spars = [3 3 -1/3 -0.2/3 0 0.2/3 -1]';

[thmm,phmm] = ode45(@lincubode,tfine,[-1 1],odeopts,spars);
phmm = {phmm(:,1) phmm(:,2)};

coefs = getcellcoefs(smooths);
newcoefs = lsqnonlin(@SplineCoefErr,coefs,[],[],...
    lsopts_out,basis_cell,Ycell,Tcell,wts,lambda,fn,[],spars);

%% Perform the Profiled Estimation
%
% |Profile_GausNewt| runs the Guass-Newton iteration for the outer
% optimization in profiling. It outputs the new parameter estimates along
% with a cell-array of functional data objects that give the model-based
% smooth to the data 

[newpars,newDEfd2] = Profile_GausNewt(startpars,lsopts_out,DEfd,fn,...
    lambda,Ycell,Tcell,wts,[],lsopts_in);

disp(['New parameter estimates: ',num2str(newpars')]);

% plot smooth with profile-estimated parameters

figure(4)
devals = eval_fdcell(tfine,newDEfd2,0);
for i = 1:length(Ycell)
    subplot(length(Ycell),1,i)
    plot(tfine,devals{i},'r','LineWidth',2);
    hold on;
%    plot(Tcell{i},Ycell{i},'b.');
    plot(plot_time,plot_path(:,i),'g','linewidth',2);
    hold off;
    if i==1
        ylabel('\fontsize{13} V')
        title(['\fontsize{13} Raw data (.), ', ...
               'profiled solution (r-) and true path (g-)'])
    else
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} R')
    end
end



function Fpermstr = Fperm_fd(yfdPar, xfdcell, betacell, wt, ...
                             nperm, argvals, q, plotres,qvec)
%  FPERM_FD does a permutation test for a functional parameter effect.
%  Arguments;
%  The first four are same as the corresponding arguments for function
%  fRegress; yfdpar, xfdCell, betacell, wt
%  NPERM   ... number of permutations to be used
%  ARGVALS ... argument values for evaluating functional responses,
%  Q       ... tail probability of quantile to compare critical value pointwise statistic
%  PLOTRES ... If nonzero, results are plotted
%  Qvec       ... tail probability of quantile to compare critical value for the maximum of the statistic
%  Results;
%  FPERMSTR;
%  A struct object with slots;
%
% Javier Apfeld  changes:
% added multiple lines for a set of q values

%  Last modified 22 September 2009 by Jim Ramsay



%  Set default arguments

if nargin < 8, plotres = 1;  end
if nargin < 7, q = 0.05;     end

if nargin < 9, qvec = q;  end %JA

if nargin < 6, argvals = []; end
if nargin < 5, nperm = 200;  end
if nargin < 4, wt = [];      end

%

q = 1-q;
qvec = 1-qvec;

if isa_fd(xfdcell{1})
    N = size(getcoef(xfdcell{1}),2);
else
    N = length(xfdcell{1});
end

tic;
fRegressStr = fRegress(yfdPar, xfdcell, betacell, wt);
elapsed_time = toc;

if elapsed_time > 30/nperm
    disp(['Estimated Computing time = ', ...
        num2str(nperm*elapsed_time),' seconds.'])
end

yhatfd = fRegressStr.yhat;

[Fvals, argvals] = Fstat_fd(yfdPar,yhatfd,argvals);
nargs = length(argvals)

Fobs = max(Fvals)

Fnull     = zeros(nperm,1);

Fnullvals = [];

for i = 1:nperm
    yfdPari      = yfdPar(randperm(N));
    fRegressStri = fRegress(yfdPari, xfdcell, betacell, wt);
    yhatfdi      = fRegressStri.yhat;
    Fi           = Fstat_fd(yfdPar,yhatfdi,argvals);
    Fnullvals    = [Fnullvals,Fi];
    Fnull(i) = max(Fnullvals(:,i));
end

pval = mean(Fobs < Fnull);
qval = quantile(Fnull,q);
qval2 = quantile(Fnull,qvec);
Fnullvals= Fnullvals'%JA
%size(Fnullvals)%
%size(Fvals)%

pvals_pts = zeros(nargs,1);
qvals_pts = pvals_pts;
for i=1:nargs
    pvals_pts(i) = mean(Fvals(i) < Fnullvals(:,i));
    qvals_pts(i) = quantile(Fnullvals(:,i),q);
end

legend_names =[];
color_list = cptcmap('CM_Paired_08','ncol',length(qvec));
if plotres
    if isa_fd(yfdPar)
        plot(argvals, Fvals, 'b-', ...
            argvals, qvals_pts, 'b--')
        legend_names={'Observed Statistic', ['pointwise ',num2str(1-q),' critical value']};
        
        for i = 1:length(qvec)
            line([min(argvals),max(argvals)],[qval2(i) qval2(i)],'Color', color_list(i,:),'LineStyle',':');
            legend_names{i+2} = ['maximum ',  num2str(1-qvec(i)),' critical value'];
        end
        xlabel('\fontsize{13} argument values')
        ylabel('\fontsize{13} F-statistic')
        title('\fontsize{13} Permutation F-Test')
        legend(legend_names)

    else
        cnts  = hist(Fnull);
        xvals = [Fnull;Fobs];
        xmax  = ceil(max(xvals));
        xmin  = min(xvals);
        Nmax  = max(cnts);
        hist(Fnull)
        hold on
        phdl=plot([Fobs, Fobs], [0,Nmax], 'b-', ...
                  [qval, qval], [0,Nmax], 'b--');
        set(phdl, 'LineWidth', 2)

        hold off
        axis([xmin,xmax,0, max(N)])
        xlabel('\fontsize{13} F-value')
        title('\fontsize{13} Permutation F-Test')
        legend('Frequency', 'Observed Statistic', ...
            ['Permutation ',num2str(1-q),' critical value'],  'location', 'EastOutside')  )        
    end
end

 

Fpermstr.pval        = pval;
Fpermstr.qval        = qval;
Fpermstr.qval2       = qval2;
Fpermstr.Fobs        = Fobs;
Fpermstr.Fnull       = Fnull;
Fpermstr.Fvals       = Fvals;
Fpermstr.Fnullvals   = Fnullvals;
Fpermstr.pvals_pts   = pvals_pts;
Fpermstr.qvals_pts   = qvals_pts;
Fpermstr.fRegressStr = fRegressStr;
Fpermstr.argvals     = argvals;


function tpermStr = tperm_fd(x1fd, x2fd, nperm, q, argvals, plotres,qvec)
%  TPERM computes a permutation t-test of the difference between
%  the functional means of two groups of data.
%
%  Arguments:
%  X1FD    ... First  sample of functional observations
%  X2FD    ... Second sample of functional observations
%  NPERM   ... Number of permutations ... default 200
%  Q       ... Tail probability for test ... default 0.05
%  ARGVALS ... Argument values at which to evaluate functions
%  PLOTRES ... If nonzero, plot results ... default 1
%  Qvec    ... tail probability of quantile to compare critical value for the maximum of the statistic
%  Results:
%
%  Javier Apfeld  changes:
%  added multiple lines for a set of q values

%  Last modified 10 October 2009 by Jim Ramsay

%  check first two arguments

if nargin < 2
    error('Arguments X1FD and X2FD not supplied.');
end

if ~isa_fd(x1fd) || ~isa_fd(x2fd)
    error('x1fd and x2fd must both be functional data objects');
end

if length(size(getcoef(x1fd))) > 2 || length(size(getcoef(x2fd))) > 2
    error('Both of X1FD and X2FD are not univariate.');
end

%  Set default arguments

if nargin < 6,  plotres = 1;   end
if nargin < 5,  argvals = [];  end
if nargin < 4,  q = 0.05;      end

if nargin < 7, qvec = q;  end %JA

if nargin < 3,  nperm = 200;   end


range1 = getbasisrange(getbasis(x1fd));
range2 = getbasisrange(getbasis(x2fd));


if ~all(range1 == range2)
    error('x1fd and x2fd do not have the same range.');
end

if isempty(argvals)
    narg = 101;
    argvals = linspace(range1(1),range1(2),narg)';
else
    narg = length(argvals);
end

q = 1-q;
qvec = 1-qvec;

x1mat = eval_fd(argvals,x1fd);
x2mat = eval_fd(argvals,x2fd);

Xmat = [x1mat,x2mat];

n1 = size(x1mat,2);
n2 = size(x2mat,2);

Tnull = zeros(nperm,1);

Tnullvals = zeros(length(argvals),nperm);

for i = 1:nperm
    tXmat = Xmat(:,randperm(n1+n2));
    
    tmean1 = mean(tXmat(:,1:n1),     2);
    tmean2 = mean(tXmat(:,n1+(1:n2)),2);
    
    tvar1 = var(tXmat(:,1:n1),     0,2)/n1;
    tvar2 = var(tXmat(:,n1+(1:n2)),0,2)/n2;
    
    Tnullvals(:,i) = abs(tmean1-tmean2)./sqrt(tvar1+tvar2);
    Tnull(i) = max(Tnullvals(:,i));
end

mean1 = mean(Xmat(:,1:n1),     2);
mean2 = mean(Xmat(:,n1+(1:n2)),2);

var1 = var(Xmat(:,1:n1),     0,2)/n1;
var2 = var(Xmat(:,n1+(1:n2)),0,2)/n2;

Tvals = abs(mean1-mean2)./sqrt(var1+var2);
Tobs  = max(Tvals);

pval = mean( Tobs < Tnull );
qval = quantile(Tnull, q);
qval2 = quantile(Tnull,qvec);

pvals_pts = zeros(narg,1);
qvals_pts = pvals_pts;
for i=1:narg
    pvals_pts(i) = mean(Tvals(i) < Tnullvals(i,:));
    qvals_pts(i) = quantile(Tnullvals(i,:), q);
end

legend_names =[];
color_list = cptcmap('CM_Paired_08','ncol',length(qvec));

if plotres
    
    ylim = [ min([Tvals;qvals_pts]),max([Tobs;qval])];
    
    plot(argvals, Tvals, 'b-', argvals, qvals_pts, 'b--')
    legend_names={'Observed Statistic', ['pointwise ',num2str(1-q),' critical value']};
    hold on
      for i = 1:length(qvec)
            line([min(argvals),max(argvals)],[qval2(i) qval2(i)],'Color', color_list(i,:),'LineStyle',':');
            legend_names{i+2} = ['maximum ',  num2str(1-qvec(i)),' critical value'];
        end
    hold off
      
    xlabel('\fontsize{13} argument values')
    ylabel('\fontsize{13} t-statistic')
    axis([min(argvals),max(argvals),ylim])
    legend(legend(legend_names), ...
           'location', 'EastOutside')  
end

tpermStr.pval         = pval;
tpermStr.qval         = qval;
tpermStr.qval2        = qval2;
tpermStr.Tobs         = Tobs;
tpermStr.Tnull        = Tnull;
tpermStr.Tvals        = Tvals;
tpermStr.Tnullvals    = Tnullvals;
tpermStr.pvals_pts    = pvals_pts;
tpermStr.qvals_pts    = qvals_pts;
tpermStr.argvals      = argvals;


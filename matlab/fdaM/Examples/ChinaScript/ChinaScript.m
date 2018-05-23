addpath ('c:\Program Files\matlab\fdaM\examples\ChinaScript')

%  load the data from file ChinaScript.mat

load ChinaScript

penpos = ChinaScript.penpos;

mark = 1:12:601;

i    = 1;
StatSci1 = penpos(:,i,:);
% Where does the pen leave the paper?
thresh = quantile(StatSci1(:, 3), .8);

sel1 = (StatSci1(:,3) < thresh);
StatSci1(StatSci1(:,3) < thresh,1:2) = NaN;
plot(StatSci1(:,1), StatSci1(:,2), 'b-', ...
     StatSci1(mark, 1), StatSci1(mark, 2), 'o')
xlabel('\fontsize{13} X') 
ylabel('\fontsize{13} Y')

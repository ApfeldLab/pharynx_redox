%
% Ramsay, Hooker & Graves (2009)
% Functional Data Analysis with R and Matlab (Springer)
%

%  Remarks and disclaimers

%  These R commands are either those in this book, or designed to 
%  otherwise illustrate how R can be used in the analysis of functional
%  data.  
%  We do not claim to reproduce the results in the book exactly by these 
%  commands for various reasons, including:
%    -- the analyses used to produce the book may not have been
%       entirely correct, possibly due to coding and accuracy issues
%       in the functions themselves 
%    -- we may have changed our minds about how these analyses should be 
%       done since, and we want to suggest better ways
%    -- the R language changes with each release of the base system, and
%       certainly the functional data analysis functions change as well
%    -- we might choose to offer new analyses from time to time by 
%       augmenting those in the book
%    -- many illustrations in the book were produced using Matlab, which
%       inevitably can imply slightly different results and graphical
%       displays
%    -- we may have changed our minds about variable names.  For example,
%       we now prefer "yearRng" to "dayrange" for the weather data.
%    -- three of us wrote the book, and the person preparing these scripts
%       might not be the person who wrote the text
%  Moreover, we expect to augment and modify these command scripts from time
%  to time as we get new data illustrating new things, add functionality
%  to the package, or just for fun.

%
% ch. 11  Functional Models and Dynamics
%

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = '../Matlab/fdaM';
addpath(fdaMPath)

%  Path to the folder containing the examples

examplesPath = [fdaMPath,'/examples'];
addpath(examplesPath)

%
% Section 11.1  Introduction to Dynamics
%
%  Figure 11.1
%

%  Purely schematic, and not reproduced in Matlab.

%
% Section 11.2 Principal Differential Analysis for Linear Dynamics
%

%  (no computations in this section)

%
% Section 11.3 Principal Differential Analysis of the Lip Data
%

lipPath = [examplesPath,'/lip'];
addpath(lipPath)

load lip

% 1.  Create an 'fd' object 'lipfd'
%
%  ----------  set up the b-spline basis object  ------------
%       use order 6 splines so we can look at acceleration

lipbasis = create_bspline_basis([0,0.35], 31, 6);

%  -----  create the fd object  ---------

%  set up the functional parameter object

Lfdobj = int2Lfd(4);  %  penalize fourth derivative
lambda = 1e-12;       %  use light smoothing
lipfdPar = fdPar(lipbasis, Lfdobj, lambda);

%  carry out the smoothing
 
lipfd = smooth_basis(liptime, lip, lipfdPar);

%  add names to dimensions

lipfd_fdnames{1} = 'Normalized time (seconds)';
lipfd_fdnames{2} = 'Replications';
lipfd_fdnames{3} = 'Lower lip position (mm)';
lipfd = putnames(lipfd, lipfd_fdnames);

% Figure 11.2

subplot(1,1,1)
plot(lipfd)

%  pda of lip data

bwtcell = cell(1,1,2);
bwtcell{1,1,1} = fdPar(lipbasis,2,0);
bwtcell{1,1,2} = fdPar(lipbasis,2,0);

xfdcell{1} = lipfd;

lipbwtestcell = pda_fd(xfdcell, bwtcell);

dfd = 0.25.*getfd(lipbwtestcell{2})^2 - getfd(lipbwtestcell{1});
             
fdnames = {'time','rep','discriminant'};

dfd = putnames(dfd, fdnames);

% Figure 11.3

subplot(3,1,1)
bfd0 = getfd(lipbwtestcell{1});
plot(liptime, eval_fd(liptime, bfd0), '-', [0,0.35], [0,0], 'r:')
ylabel('beta 0')
subplot(3,1,2)
bfd1 = getfd(lipbwtestcell{2});
plot(liptime, eval_fd(liptime, bfd1), '-', [0,0.35], [0,0], 'r:')
ylabel('beta 1')
subplot(3,1,3)
plot(liptime, eval_fd(liptime, dfd), '-', [0,0.35], [0,0], 'r:')
xlabel('time')
ylabel('discriminant')

% Figure 11.4

pda_overlay(bwtcell)

%
% Section 11.4 PDA of the Handwriting Data
%

handwritePath = [examplesPath,'/handwrit'];
addpath(handwritePath)

fid      = fopen('fdareg.dat','rt');
fdaarray = reshape(fscanf(fid,'%f'), [20,2,1401]);
fdaarray = permute(fdaarray,[3,1,2]);
fdaarray = fdaarray/1000;   %  convert spatial unit to meters

%  Set up time values and range.
%  It is best to choose milliseconds as a time scale
%  in order to make the ratio of the time
%  unit to the inter-knot interval not too
%  far from one.  Otherwise, smoothing parameter values
%  may be extremely small or extremely large.

fdatime  = linspace(0, 2300, 1401)';
fdarange = [0, 2300];

nbasis   = 1406;
norder   =    7;
fdabasis = create_bspline_basis(fdarange, nbasis, norder);

fdafd  = fd(zeros(nbasis,20,2), fdabasis);
lambda = 1e9;
fdaPar = fdPar(fdafd, 5, lambda);

%  set up the functional data structure 

[fdafd, df, gcv] = smooth_basis(fdatime, fdaarray, fdaPar);
%  Add suitable names for the dimensions of the data.
fdafd_fdnames{1} = 'Milliseconds';
fdafd_fdnames{2} = 'Replications';
fdafd_fdnames{3} = 'Metres';
fdafd = putnames(fdafd, fdafd_fdnames);

xfdcell = cell(2,1);
xfdcell{1} = fdafd(:,1);
xfdcell{2} = fdafd(:,2);

pdaPar = fdPar(fdabasis,2,1);

bwtcell = cell(2,2,2);
for i=1:2
    for j=1:2
        for k=1:2
            bwtcell{i,j,k} = pdaPar;
        end
    end
end
          
bwtcell = pda_fd(xfdcell, bwtcell);

% Figure 11.5

eigen_pda(bwtcell);

%
% Section 11.5 Registration and PDA
%

WfdPar = fdPar(getbasis(lipfd),2,1e-16);
lipmarks = lipmarks*0.35;
lipmeanmarks= mean(lipmarks);

[lipregfd, warpfd, Wfd] = ...
    landmarkreg(lipfd, lipmarks, lipmeanmarks, WfdPar);
D1lipfd = deriv_fd(lipfd,1);
D2lipfd = deriv_fd(lipfd,2);
D1lipregfd = register_fd(mean(D1lipfd), D1lipfd, warpfd);
D2lipregfd = register_fd(mean(D2lipfd), D2lipfd, warpfd);
xfdcell{1} = -D1lipregfd;
xfdcell{2} = -lipregfd;

bwtcell = cell(1,2);
bwtcell{1,1} = fdPar(lipbasis,2,0);
bwtcell{1,2} = fdPar(lipbasis,2,0);
xfdcell = cell(1,1);
xfdcell{1} = lipfd;
lipbwtestcell1 = pda_fd(xfdcell, bwtcell);

lipregpdaStr = fRegress(D2lipregfd, xfdcell, bwtcell);

lipbwtestcell2 = lipregpdaStr.betahat;

% Figure 11.6

subplot(2,1,1)
bwt11fd = getfd(lipbwtestcell2{1});
phdl1 = plot(liptime, eval_fd(liptime, bwt1fd));
set(phdl1,'LineWidth', 2)
hold on
bwt12fd = getfd(lipbwtestcell1{1});
phdl2 = plot(liptime, eval_fd(liptime, bwt12fd));
set(phdl2,'LineWidth', 2, 'LineStyle', '--')
hold off
subplot(2,1,2)
bwt21fd = getfd(lipbwtestcell2{2});
phdl1 = plot(liptime, eval_fd(liptime, bwt21fd));
set(phdl1,'LineWidth', 2)
hold on
bwt22fd = getfd(lipbwtestcell1{2});
phdl2 = plot(liptime, eval_fd(liptime, bwt22fd));
set(phdl2,'LineWidth', 2, 'LineStyle', '--')
hold off

%
% Section 11.6 Details for pda_fd, eigen_fd, pda_overlay
%              and register_newfd
%



%
% Section 11.7 Some Things to Try
%
% (exercises for the reader)

%
% Section 11.8  More to Read
%

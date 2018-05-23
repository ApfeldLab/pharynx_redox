%%%   Smoothing Operators in the FDA Package
%
%%   Last Modified: 14/02/2007
%
%%   This module introduces the FDA functionality for smoothing data using
%%   linear differential operator penalties. 
%
%%%   Basis Setup
%
%%   We will smooth using the handwriting data

load 'vancprec.mat'
%%

daytime = vancprec(:,1);
prec = vancprec(:,2);

plot(daytime,prec,'.')

%%   Now set up a basis

basis_range = [min(daytime) max(daytime)];

norder = 5;

breaks = linspace(basis_range(1),basis_range(2),300);

nbasis = length(breaks) + norder - 2;

bspline_obj = create_bspline_basis(basis_range,nbasis,norder,breaks);
%%

plot(bspline_obj)
%

%%   And we might even set this to data

precfd = data2fd(daytime,prec,bspline_obj);
%%

prechat = eval_fd(daytime,precfd);
figure
plot(daytime,prechat)
hold on
plot(daytime,prec,'.')
hold off

% %%   Defining Linear Differential Operators
% 
% %   In order to define a penalty, a linear differential operator is
% %   necessary. The FDA package has a class LFD to do this.
% 
% %   Formally, an linear differential operator takes the form
% 
% %%%      Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) + \exp[w_m(t)]
% %%%      D^m x(t)
% 
% %   and we can define the w_i(t) as functional data objects. 
% 
% %   To create an Lfd Object, we call a function Lfd. this takes two arguments
% 
% %   1. m - the highest derivative in L
% 
% %   2. wfdcell - a cell array containing the w functions. 
% 
% %   See 

help Lfd

%%   Setting up all of that is rather annoying, so there are some shortcuts. 
% 
% %%   Smoothing by derivatives
% 
% %   If I only want to penalize D^m then I can use int2Lfd.
% 
% %   For the standard second derivative penalty:

Lfd2 = int2Lfd(2)

%%   Functional Parameter Objects
% 
% %   And fdPar object contains a lot of objects and avoids needing to type
% %   them all as arguments to functions. I'm not convinced that this is
% %   always helpfull, but the code is set up that way.

lambda =200;

precPar = fdPar(bspline_obj,Lfd2,lambda);

precfd = smooth_basis(daytime,prec,precPar);

precvals = eval_fd(daytime,precfd);

plot(daytime,precvals)
hold on
plot(daytime,prechat,'r--')
%plot(daytime,prec,'.')
hold off
%%
figure
plot(daytime,prechat-precvals)

%%%   GCV and DF
% 
% %   I can get both the equivalent degrees of freedom and GCV score from this,
% %   too

[precfd,df,gcv] = smooth_basis(daytime,prec,precPar);

df
gcv

%%%   Trying out the Weather Data
%
%%   Remember, typically, lambda varies on the log scale. 

df = zeros(12,1);
gcv = zeros(12,1);

for i = -1:10
    
    j = i + 2;
    subplot(2,6,j)
    lambda = 10^i;

    precPar = fdPar(bspline_obj,Lfd2,lambda);  

    [precfd,df(j),gcv(j)] = smooth_basis(daytime,prec,precPar);

    precvals = eval_fd(daytime,precfd);

    plot(daytime,precvals,'linewidth',2)
    hold on
   % plot(daytime,prec,'r.')
    hold off
    title(strcat('log lambda = ',int2str(i),'gcv=',num2str(gcv(j))),'fontsize',10);
   

end

%%   Now I can look at

df

plot(-1:10,gcv,'r')

[m,i] = min(gcv)

%%   So choose the minimizing value

lambda = 10^3;

precPar = fdPar(bspline_obj,Lfd2,lambda);

[precfd,df(j),gcv(j)] = smooth_basis(daytime,prec,precPar);

precvals = eval_fd(daytime,precfd);

plot(daytime,precvals,'linewidth',2)
hold on
plot(daytime,prec,'r.')
hold off
%%

% %   Note that smooth_basis also accepts a vector of weights as a fourth
% %   argument and a multiplication factor as a fifth argument. 
% 
% %   Exercize - try the above with a multiplication factor of 1.4. 
% 
% %   See also smooth_fd. 
% 
% %%   More Complex Penalties
% 
% %   Constant co-efficient differential equations are the next most common
% %   penalty. In particular, the harmonic acceleration penalty is frequently
% %   used for data with a known period. This penalty is 
% 
% %%%                   D^3 + a^2 D
% 
% %   In this case our co-efficient functions are constant with values 
% %   [0 a^2 0 1] but since the last one is exponentiated, we leave it off. 

Lbasis  = create_constant_basis([0,365]);  %%%%      create a constant basis
Lcoef   = [0,(2*pi/365)^2,0];    %%%%      set up three coefficients
wfd     = fd(Lcoef,Lbasis);      %%   define an FD object for weight functions
wfdcell = fd2cell(wfd);          %%   convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %%%%      define the operator object

%%   In particular, fd2cell turns a functional data object with replications
%   into a cell array that can then be used in harmaccelLfd.

wfdcell

%%   Now we can try the same thing


figure;
df = zeros(7,1);
gcv = zeros(7,1);

for i = -1:5
    
    j = i + 2;
    
    lambda = 10^(2*i);

    precPar = fdPar(bspline_obj,harmaccelLfd,lambda);

    [precfd,df(j),gcv(j)] = smooth_basis(daytime,prec,precPar);

    precvals = eval_fd(daytime,precfd);

    plot(daytime,precvals,'linewidth',2)
    hold on
    plot(daytime,prec,'r.')
    hold off
    title(strcat('log lambda = ',int2str(2*i)),'fontsize',20);
    pause

end

%%   Now I can look at

df

plot(-1:2:11,gcv,'r')

[m,i] = min(gcv)

%%   So choose the minimizing value

lambda = 10^6;

precPar = fdPar(bspline_obj,harmaccelLfd,lambda);

[precfd,df(j),gcv(j)] = smooth_basis(daytime,prec,precPar);

precvals = eval_fd(daytime,precfd);

plot(daytime,precvals,'linewidth',2)
hold on
plot(daytime,prec,'r.')
hold off

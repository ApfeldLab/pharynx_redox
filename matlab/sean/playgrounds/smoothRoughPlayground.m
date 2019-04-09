i = 41;

smoothLambda = 10^0.0891;
roughLambda  = 10^-3;

smoothF1 = makeWormFd_SJ(unreg_i1(:,i), 'lambda', smoothLambda);
roughF1 = makeWormFd_SJ(unreg_i1(:,i), 'lambda', roughLambda);

smoothF2 = makeWormFd_SJ(unreg_i2(:,i), 'lambda', smoothLambda);
roughF2 = makeWormFd_SJ(unreg_i2(:,i), 'lambda', roughLambda);

figure;
hold on;
scatter(1:100, unreg_i1(:,i));
plot(smoothF1);
plot(roughF1);
legend('Data', 'Smooth', 'Rough');

Wfd0 = fdObjs(i).warpFD(2);

%  set up a fine mesh of argument values
ybasis  = getbasis(smoothF1);
ynbasis = getnbasis(ybasis);
nfine = max([201,10*ynbasis + 1]);

wcoef  = getcoef(Wfd0);
wbasis = getbasis(Wfd0);
wtype  = getbasistype(wbasis);
nbasis = getnbasis(wbasis);
norder = nbasis - length(getbasispar(wbasis));
rangex = getbasisrange(wbasis);
wdim   = size(wcoef);
ncoef  = wdim(1);

xlo   = rangex(1);
xhi   = rangex(2);
width = xhi - xlo;
xfine = linspace(xlo, xhi, nfine)';

penmat  = eval_penalty(ybasis);
penmat  = penmat + 1e-10 .* max(max(penmat)) .* eye(ynbasis);
penmat  = sparse(penmat);

yfine = squeeze(eval_fd(xfine, roughF2));


% JMAX = 15;
% basiscell = cell(1,JMAX);
% 
% ffine    =   monfn(xfine, Wfd0, basiscell);
% fmax     = ffine(nfine);
% hfine    = xlo + width.*ffine./fmax;
hfine = eval_fd(xfine, Wfd0);
hfine(1)     = xlo;
hfine(nfine) = xhi;

roughRegF2 = regyfn(xfine, yfine, hfine, roughF2, Wfd0, penmat, 0);

figure;
hold on;
scatter(1:100, unreg_i1(:,i));
scatter(1:100, unreg_i2(:,i));
plot(roughF1);
plot(roughF2);
plot(roughRegF2);
legend('Data1', 'Data2', 'Rough F1', 'Rough F2', 'Rough Reg F2');


function yregfd = regyfn(xfine, yfine, hfine, yfd, Wfd, penmat, ...
                         periodic)

%  get shift value for the periodic case from Wfd

coef   = getcoef(Wfd);
shift  = coef(1);  
coef(1)= 0;

%  if all coefficients are zero, no transformation  needed

if all(coef == 0)
   if periodic
      if shift == 0
         yregfd = yfd;
         return;
      end
   else
      yregfd = yfd;
      return;
   end
end

%  Estimate inverse of warping function at fine mesh of values  
%  28 dec 000
%  It makes no real difference which 
%     interpolation method is used here.
%  Linear is faster and sure to be monotone.
%  Using WARPSMTH added nothing useful, and was abandoned.

nfine       = length(xfine);
hinv        = safeinterp(hfine, xfine, xfine);
hinv(1)     = xfine(1);
hinv(nfine) = xfine(nfine);

%  carry out shift if period and shift ~= 0

if periodic && shift ~= 0
   yfine = shifty(xfine, yfine, shift);
end

%  smooth relation between Y and HINV
%  this is the same code as in PROJECT_BASIS, but avoids
%  recomputing the penalty matrix

basis    = getbasis(yfd);
if any(hinv < xfine(1))
    save hfine
    save xfine
    save hinv
    error(['HINV values too small by ',num2str(xfine(1)-min(hinv))]);
end
if any(hinv > xfine(nfine))
    error(['HINV values too large by ',num2str(max(hinv)-xfine(nfine))]);
end
basismat = getbasismatrix(hinv, basis);
Bmat     = basismat' * basismat;
lambda1  = (0.0001 .* sum(diag(Bmat)))./sum(diag(penmat));
Cmat     = Bmat + lambda1 .* penmat;
Dmat     = basismat' * yfine;
ycoef    = symsolve(Cmat,Dmat);

%  set up FD object for registered function

yregfd   = fd(ycoef, basis);
end

%%
figure;
hold on;
scatter(1:100, unreg_i1(:,i));
plot(fdObjs.);
plot(roughF1);
legend('Data', 'Smooth', 'Rough');

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
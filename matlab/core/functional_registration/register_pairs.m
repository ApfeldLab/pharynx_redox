function [reg_fd_A, reg_fd_B] = register_pairs(A, B, ...
    f_n_basis, f_order, f_deriv_penalty, f_lambda,...
    warp_n_basis, warp_order, warp_deriv_penalty, warp_lambda)

A = A';
B = B';

f_n_basis = 64;
f_order = 4;
f_deriv_penalty = 2;
f_lambda = 10^-.5;

xs = linspace(0, 100, size(A, 1));
f_basis = create_bspline_basis([0, 100], f_n_basis , f_order);
f_fd_par = fdPar(f_basis, f_deriv_penalty, f_lambda);

fd_A_der = deriv_fd(smooth_basis(xs, A, f_fd_par), 1);
fd_B_der = deriv_fd(smooth_basis(xs, B, f_fd_par), 1);

fd_A = smooth_basis(xs, A, f_fd_par);
fd_B = smooth_basis(xs, B, f_fd_par);

warp_n_basis = 3;
warp_order = 3;
warp_deriv_penalty = int2Lfd(2);
warp_lambda = 5e3;

warp_basis = create_bspline_basis([0, 100], warp_n_basis, warp_order);
warp_fd_par = fdPar(warp_basis, warp_deriv_penalty, warp_lambda);

[~, warp_fd, wfd] = register_fd(fd_A_der, fd_B_der, warp_fd_par);

nfine = 1000;
xfine = linspace(1,100,nfine)';
yfine = eval_fd(xfine, fd_B);

ybasis  = getbasis(fd_B_der);
ynbasis = getnbasis(ybasis);
ybasismat = cell(1,4);
ybasismat{1,1} = xfine;
for ideriv=0:2
    ybasismat{1,ideriv+2} = eval_basis(xfine, ybasis, ideriv);
end
ybasis = putbasisvalues(ybasis, ybasismat);
penmat  = eval_penalty(ybasis);
penmat  = penmat + 1e-10 .* max(max(penmat)) .* eye(ynbasis);
penmat  = sparse(penmat);

basiscell = cell(1,15);
ffine    =   monfn(xfine, warp_fd, basiscell);
fmax     = ffine(nfine);
width = 100 - 1;
hfine    = 1 + width.*ffine./fmax;
hfine(1)     = 1;
hfine(nfine) = 100;

reg_fd_A = fd_A;
for i=1:size(A, 1)
    reg_fd_B(i) = regyfn(xfine, yfine(:,i), hfine(:,i), fd_B(i), wfd, penmat, 0);
end
% reg_fd_A = eval_fd(xs, fd_A);
% reg_fd_B = eval_fd(xs, reg_fd_B);
end
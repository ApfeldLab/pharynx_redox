function [reg_fd_A, reg_fd_B] = register_pairs(A, B, ...
    f_n_basis, f_order, f_deriv_penalty, f_lambda,...
    warp_n_basis, warp_order, warp_deriv_penalty, warp_lambda)

A = A';
B = B';

% n_basis = 64;
% order = 4;
% deriv_penalty = 2;
% lamda = 1e-6;

xs = linspace(0, 100, size(A, 1));
f_basis = create_bspline_basis([0, 100], f_n_basis , f_order);
f_fd_par = fdPar(f_basis, f_deriv_penalty, f_lambda);

fd_A = smooth_basis(xs, A, f_fd_par);
fd_B = smooth_basis(xs, B, f_fd_par);

% warp_n_basis = 6;
% warp_order = 6;
% warp_deriv_penalty = int2Lfd(4);
% warp_lambda = 5e3;

warp_basis = create_bspline_basis([0, 100], warp_n_basis, warp_order);
warp_fd_par = fdPar(warp_basis, warp_deriv_penalty, warp_lambda);

[reg_fd_B, ~, ~] = register_fd(fd_A, fd_B, warp_fd_par);

reg_fd_A = eval_fd(xs, fd_A);
reg_fd_B = eval_fd(xs, reg_fd_B);
end
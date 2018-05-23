% Use GCV Minimization Criteria to find optimal Lambda (smoothing
% parameter)

%% Load Data
p_410_m_410 = csvread(...
    '../experiments/single_v_dual_poly/poly_410_measure_410_sean.csv', ...
    1, 1);
p_410_m_410(p_410_m_410<500) = 0;

sq = flipud(ssquare(p_410_m_410));

subset = 1:10;
data = sq(:,subset);
positions = 1:size(sq,1);

initialFd = makeWormFd_SJ(data, 'lambda', 10^-2);

%% Lambda Optimization
log_lambdas = -10:0.01:1;

Lfd2 = int2Lfd(2);

gcv_vals = zeros(length(log_lambdas),1);
df_vals = gcv_vals;
tic
parfor i=1:length(log_lambdas)
    lambda_i = 10.^log_lambdas(i);
    fdParams_i = fdPar(getbasis(initialFd), Lfd2, lambda_i);
    
    [fd_i, dftmp, gcvtmp_i] = smooth_basis(positions, data, fdParams_i);
    gcv_vals(i) = sum(gcvtmp_i);
    df_vals(i) = dftmp;
    vals_i = eval_fd(positions, fd_i);
end
toc

[min_gcv, min_gcv_idx] = min(gcv_vals);
min_loglambda = log_lambdas(min_gcv_idx);

figure('Name', 'GCV vs log lambda');
plot(log_lambdas, gcv_vals); hold on;
scatter(min_loglambda, min_gcv);
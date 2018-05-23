function [optimalLambda, gcv_vals] = findOptimalLambda(sqData, varargin)
% findOptimalLambda Find the optimal Smoothing Parameter (lambda) for the 
% given data.
%   REQUIRED:
%       sqData:
%           This should be an MxN matrix (M=length of intensity vector, 
%           N=number of animals).
%   OPTIONAL:
%       logLambdas:
%           1xN matrix, each column consisting of the log(lambda) to test.
%           Can create with: minLambda:step:maxLambda or linspace(minLambda, maxLambda, N)
%           DEFAULT: -10:0.01:1
%       fdObj:
%           An fd object, made from sqData. Overrides the `basis` parameter
%       basis:
%           A fdBasis object, capable of representing sqData.
%           DEFAULT: bspline basis with n_basis=96, n_order=6
%       plotGCV:
%           0/1, indicating whether or not to plot GCV/log(lambda)

    LOG_LAMBDA_DEFAULT = -10:0.01:1;
    PLOT_DEFAULT = 0;
    BASIS_DEFAULT = create_bspline_basis([1 size(sqData, 1)], 96, 6);
    DIFF_PENALTY_DEFAULT = 2;
 
    persistent p;
    if isempty(p)
        p = inputParser;
        p.FunctionName = 'findOptimalLambda';
        
        addParameter(p, 'logLambdas', LOG_LAMBDA_DEFAULT);
        addParameter(p, 'basis', BASIS_DEFAULT);
        addParameter(p, 'plotGCV', PLOT_DEFAULT);
        addParameter(p, 'diffPenalty', DIFF_PENALTY_DEFAULT);
    end
    parse(p,varargin{:});
    
    log_lambdas = p.Results.logLambdas;
    basis_ = p.Results.basis;
    positions = 1:size(sqData, 1);

    Lfd2 = int2Lfd(p.Results.diffPenalty);

    gcv_vals = zeros(length(log_lambdas),1);

    parfor i=1:length(log_lambdas)
        lambda_i = 10.^log_lambdas(i);
        fdParams_i = fdPar(basis_, Lfd2, lambda_i);

        [~, ~, gvc_i] = smooth_basis(positions, sqData, fdParams_i);
        gcv_vals(i) = sum(gvc_i);
        % vals_i = eval_fd(positions, fd_i);
    end

    [min_gcv, min_gcv_idx] = min(gcv_vals);
    min_loglambda = log_lambdas(min_gcv_idx);
    optimalLambda = 10^min_loglambda;

    if p.Results.plotGCV
        figure('Name', 'GCV vs log lambda');
        plot(log_lambdas, gcv_vals, 'bo-'); hold on;
        scatter(min_loglambda, min_gcv, 'r');
        xlabel('\fontsize{13} log_{10}\lambda')
        ylabel('\fontsize{13} GCV(\lambda)')
    end
end
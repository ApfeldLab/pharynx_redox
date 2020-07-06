function smoothed_data = ...
    smooth_profiles(data, resample_resolution, N_BREAKS, ORDER, LAMBDA, N_DERIV)

    [sm_fd, ~] = makeWormFd_SJ(data, 'lambda', LAMBDA, 'n_order', ORDER, 'n_breaks', N_BREAKS);

    xs = linspace(1, 100, resample_resolution);
    smoothed_data = eval_fd(xs, sm_fd, N_DERIV);
end
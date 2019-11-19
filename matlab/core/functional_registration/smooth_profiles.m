function eval_data = smooth_profiles(profile_data, lambda, n_eval_pts)
    fdata = makeWormFd_SJ(profile_data, 'lambda', lambda);
    
    eval_data = eval_fd(linspace(1,100,n_eval_pts), fdata);
end
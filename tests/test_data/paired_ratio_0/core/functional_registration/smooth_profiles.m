function eval_data = smooth_profiles(profile_data, lambda, order, n_basis, n_eval_pts)
    fdata = makeWormFd_SJ(profile_data, 'lambda', lambda, 'n_order', order, 'n_breaks', n_basis);
    
    eval_data = eval_fd(linspace(1,100,n_eval_pts), fdata);
end
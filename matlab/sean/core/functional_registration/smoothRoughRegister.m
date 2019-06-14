function [smoothFD1, smoothFD2, roughFD1, roughFD2, regRoughFD2, regSmoothFD2, warpFD] = smoothRoughRegister(i1, i2, smoothLambda, roughLambda, warpLambda)
    % First, we create the smooth FD objects
    
    smoothFD1 = makeWormFd_SJ(i1, 'lambda', smoothLambda);
    smoothFD2 = makeWormFd_SJ(i2, 'lambda', smoothLambda);
    
    warp_nBasis = 6;
    warp_order  = 4;
    
    warpBasis = create_bspline_basis([1 100], warp_nBasis, warp_order);
    warpFDParObj = fdPar(warpBasis, int2Lfd(2), warpLambda);
    
    
    % Next, we register the second channel to the first
    
    n_animals = size(i1, 2);
    roughFD1 = makeWormFd_SJ(i1, 'lambda', roughLambda);
    roughFD2 = makeWormFd_SJ(i2, 'lambda', roughLambda);
    
    % We first register the *smooth* versions of the data
    [regSmoothFD2, warpFD] = register_fd(smoothFD1, smoothFD2, ...
        warpFDParObj, 0, 2, 1e-4, 100, 0);
    
    ybasis = getbasis(smoothFD1);
    ynbasis = getnbasis(ybasis);
    nfine = max([201,10*ynbasis + 1]);
    wbasis = getbasis(warpFD);
    rangex = getbasisrange(wbasis);
    xlo   = rangex(1);
    xhi   = rangex(2);
    xfine = linspace(xlo, xhi, nfine)';
    
    penmat  = eval_penalty(ybasis);
    penmat  = penmat + 1e-10 .* max(max(penmat)) .* eye(ynbasis);
    penmat  = sparse(penmat);
    
    regCoefs = getcoef(roughFD2);
    
    parfor i=1:n_animals
        yfine = squeeze(eval_fd(xfine, roughFD2(i)));
        hfine = eval_fd(xfine, warpFD(i));
        hfine(1)     = xlo;
        hfine(nfine) = xhi;

        roughRegF2 = regyfn(xfine, yfine, hfine, roughFD2(i), warpFD(i), penmat, 0);
        regCoefs(:,i) = getcoef(roughRegF2);
    end
    
    regRoughFD2 = fd(regCoefs, getbasis(roughFD2));
end
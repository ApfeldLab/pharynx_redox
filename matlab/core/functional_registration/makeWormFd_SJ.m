function wormFd = makeWormFd_SJ(intensityData, varargin)
    %SmoothIntensity Return a functional data object containing a smoothing of the
    %intensity Data
    %   intensityData should be length-normalized (see square).
    %   shape: (obs, profile_length)

    persistent p;

    if isempty(p)
        p = inputParser;
        p.FunctionName = 'makeWormFd_SJ';
        addParameter(p, 'lambda', .01);
        addParameter(p, 'n_order', 20);
        addParameter(p, 'n_breaks', size(intensityData, 2));
    end

    parse(p, varargin{:});
    lambda = p.Results.lambda;
    n_order = p.Results.n_order;
    n_breaks = p.Results.n_breaks;

    breaks = linspace(1, 100, n_breaks);
    n_basis = length(breaks) + n_order - 2;

    basis_range = [1 100];
    bspline_basis = create_bspline_basis(basis_range, n_basis, n_order, breaks);

    Lfd2 = int2Lfd(2);
    wormFdPar = fdPar(bspline_basis, Lfd2, lambda);
    argvals = linspace(basis_range(1), basis_range(2), size(intensityData, 2));
    [wormFd, ~, ~] = smooth_basis(argvals, intensityData.', wormFdPar);
end
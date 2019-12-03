function [reg_rgh_410, reg_rgh_470, warpfd, wfd, d410, reg_d470] = ...
    channel_register(...
    meas_410, meas_470, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    N_DERIV)

% Set the Defaults
% Chosen emperically
if nargin < 12, N_DERIV = 1          ; end
if nargin < 11, ROUGH_ORDER = 4      ; end
if nargin < 10, ROUGH_N_BASIS = 200  ; end
if nargin < 9 , ROUGH_LAMBDA = 0.1   ; end
if nargin < 8 , SMOOTH_ORDER = 4     ; end
if nargin < 7 , SMOOTH_N_BASIS = 50  ; end
if nargin < 6 , SMOOTH_LAMBDA = 2    ; end
if nargin < 5 , WARP_LAMBDA = 5      ; end
if nargin < 4 , WARP_ORDER = 4       ; end
if nargin < 3 , WARP_N_BASIS = 10    ; end
if nargin < 2 
    error('Less than three arguments supplied.'); 
end

if (~isequal(size(meas_410), size(meas_470)))
    error('Measurement arrays must be the same dimensions');
end

disp('Registering 470 to 410: ');

% Set up the intensity FD objects
sm410 = makeWormFd_SJ(meas_410, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_basis', SMOOTH_N_BASIS);
sm470 = makeWormFd_SJ(meas_470, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_basis', SMOOTH_N_BASIS);

rgh410 = makeWormFd_SJ(meas_410, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_basis', ROUGH_N_BASIS);
rgh470 = makeWormFd_SJ(meas_470, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_basis', ROUGH_N_BASIS);

% Registering the derivative works emperically better sometimes
% If N_DERIV == 0, this is just the intensity profile
d410 = deriv(sm410, N_DERIV);
d470 = deriv(sm470, N_DERIV);

% % z-standardize
% xs = linspace(0, 1, 1000);
% e410 = eval_fd(xs, d410);
% e470 = eval_fd(xs, d410);
% 
% d410 = (d410 - mean(e410, 'all')) ./ std(e410, 0, 'all');
% d470 = (d470 - mean(e470, 'all')) ./ std(e470, 0, 'all');

% Set up Warp Constraints
warpBasis = create_bspline_basis([0 1], WARP_N_BASIS, WARP_ORDER);
fdParObj = fdPar(warpBasis, int2Lfd(2), WARP_LAMBDA);

% Register the d(smooth)s to each other
periodic = 0;
crit = 2;
conv = 1e-8;
iterlim = 50;
dbglev = 1;
[reg_d470, warpfd, wfd] = register_fd(d410, d470, fdParObj, periodic, crit, conv, iterlim, dbglev);

% Use the d(smooth) to warp the rough
reg_rgh_470 = synch(linspace(0,1,1000), rgh470, wfd);
reg_rgh_410 = rgh410;
end
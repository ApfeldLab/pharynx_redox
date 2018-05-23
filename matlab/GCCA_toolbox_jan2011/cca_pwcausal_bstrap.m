function ret = cca_pwcausal_bstrap(X,Nr,Nl,nlags,nBoot,nBwin,Fs,freq,pval,CORRTYPE)
% -----------------------------------------------------------------------
%   FUNCTION: cca_pwcausal_bstrap
%   PURPOSE:  spectral granger causality with
%             bootstrap confidence interval construction
%             (for single trial use Nr=1, Nl = N)
%
%   INPUT:  X           -   nvar (rows) by nobs (cols) observation matrix
%           Nr          -   # realizations
%           Nl          -   length of each realization
%           nlags       -   number of lags to include in model
%           nBoot       -   # bootstrap resamples
%           nBwin       -   # observations in each window (if Nr>>1, set
%                           nBwin=Nl, then each trial is a separate bstrap
%                           sample)
%           Fs          -   sampling rate (Hz)
%           freq        -   vector of frequencies to analyze
%           pval        -   desired significance threshold
%           CORRTYPE    -   (0=none,1=Bonferroni,2=approx FDR)
%
%   OUTPUT: GW          -   nvar * nvar * length(freq) spectral granger
%           COH         -   nvar * nvar * length(freq) coherence
%                           (power spectrum is along diagonal)
%           pp          -   power spectrum (as above, provided separately)
%           waut        -   Durbin Watson residual autocorrelation sig value
%           cons        -   MVAR model consistency (Ding et al 2000)
%           ll          -   nvar*nvar lower CI limit on GW 
%           ul          -   nvar*nvar upper CI limit on GW 
%           ci          -   post-correction confidence interval (CI)
%
%   Based on pwcausal.m 
%
%   AKS Apr 2009
% COPYRIGHT NOTICE AT BOTTOM
% -----------------------------------------------------------------------

% check inputs
if nargin<10,
    error('error in cca_pwcausal_bstrap: insufficient inputs');
end
if nBwin == Nl & Nr<10,
    disp('warning in cca_pwcausal_bstrap: perhaps too few trials for effective bootstrapping, reduce nBwin');
end
if mod(Nl,nBwin) ~= 0,
    error('error in cca_pwcausal_bstrap: trial length must be an integer multiple of bstrap win length');
end
nwin = (Nl*Nr)/nBwin;   % figure number of bootstrap windows
if ~isint(nwin),
    error('error in cca_pwcausal_bstrap: data length must be an integer multiple of bstrap win length');
end

% figure data parameters
nobs = size(X,2);
nvar = size(X,1);

% figure confidence intervals
if CORRTYPE == 1,
    pval = pval/(nvar*(nvar-1));
elseif CORRTYPE == 2, 
     m = nvar*(nvar-1);
     pval = ((m+1)/(2*m))*pval; 
end
CI_upper = pval/2;       % upper percentile of confidence interval
CI_lower = 1-CI_upper;   % lower percentile of confidence interval
disp(['setting confidence intervals to ',num2str(CI_upper),' & ',num2str(CI_lower)]);

% remove ensemble mean
m = mean(X');
mall = repmat(m',1,nobs);
X = X-mall;

% calculate spectral GC from original data (including whiteness +
% consistency checks)
[GW,COH,pp,waut,cons] = cca_pwcausal(X,Nr,Nl,nlags,Fs,freq,1);
R_original = GW;

% generate bootstrap resamples
R_bstrp = zeros(nvar,nvar,length(freq),nBoot);       % bootstrap matrix (4D)
for ii=1:nBoot,
    XX = genBoot(X,nBwin,nwin,nvar,nobs);
    [gw2,dummy,dummy] = cca_pwcausal(XX,nwin,nBwin,nlags,Fs,freq);
    R_bstrp(:,:,:,ii) = gw2;    % need a 4D structure here
    disp(['done pwcausal bootstrap trial ',num2str(ii),'/',num2str(nBoot)]);
end

% Find the 5% and 95% empirical quantile of R(1)-R_original,....,R(nBoot)-R_original,
% Q_5 and Q_95, and obtain the 90% bootstrap confidence interval for the true R as:
% (R_original - Q_95, R_original - Q_5).
R_diff = zeros(nvar,nvar,length(freq),nBoot);
for ii=1:nBoot,
    R_diff(:,:,:,ii) = R_bstrp(:,:,:,ii)-R_original;
end

% construct confidence intervals
ll = zeros(nvar,nvar,length(freq));
ul = zeros(nvar,nvar,length(freq));
for ii=1:nvar-1,
    for jj=ii+1:nvar,
        for kk=1:length(freq),
            if ii~=jj,
                ul(ii,jj,kk) = quantile(squeeze(R_diff(ii,jj,kk,:)),CI_upper);
                ll(ii,jj,kk) = quantile(squeeze(R_diff(ii,jj,kk,:)),CI_lower);
                ul(jj,ii,kk) = quantile(squeeze(R_diff(jj,ii,kk,:)),CI_upper);
                ll(jj,ii,kk) = quantile(squeeze(R_diff(jj,ii,kk,:)),CI_lower);
            end
        end
    end
end
ul = R_original-ul;
ll = R_original-ll;
PR = ll>0;              % find significant interactions

% format output
ret.GW = GW;
ret.COH = COH;
ret.pp = pp;
ret.PR = PR;
ret.ul = ul;
ret.ll = ll;
ret.waut = waut;
ret.cons = cons;


% This file is part of GCCAtoolbox.  
% It is Copyright (C) Anil Seth (2009) 
% 
% GCCAtoolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GCCAtoolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GCCAtoolbox. If not, see <http://www.gnu.org/licenses/>.









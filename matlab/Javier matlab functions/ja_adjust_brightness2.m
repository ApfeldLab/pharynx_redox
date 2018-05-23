function [ img2 ] = ja_adjust_brightness2(img,refimg,maxcutoff,mincutoff,cm, cmin, cmax)
%ja_adjust_brightness Adjusts the brighness of one image
%   This function adjusts the brighness of an image using the intensity of
%   a reference image. 
%   Steps are as follows:
%   1) The colormapped image is transformed to HSV.
%   2) The brighness is adjusted using the intensities of the reference
%      image.  Intensities from 0 to the cutoff are adjusted linearly 
%      from 0 to 1.
%
%   img       intensity image whose brightness will be adjusted
%   refimg    reference intensity image used to adjust the brighness of img
%   maxcutoff    intensity cutoff.  set the minimal intensity were brighness is 100%
%   mincutoff    intensity cutoff.  set the maximal intensity were brighness is 0%
%   cm        colormap
%   cmin      lower bound for colormap
%   cmax      upper bound for colormap
%

if nargin < 7
    cmax = 1;
end

if nargin < 6
    cmin = 0;
end

if nargin < 5
    cm = colormap('gray');
end

if nargin < 1
    disp('missing image');
end

if nargin < 2
    refimg = ones(size(img));
end

if nargin < 3
    maxcutoff = max(max(refimg));
end

if nargin < 4
    mincutoff = min(min(refimg));
end

img_cm =[];
img_cmidx = [];
cm_length = size(cm,1);
tmp =[];

img(isnan(img)) = 0;
img(isinf(img)) = 0;

refimg(isnan(refimg)) = 0;
refimg(isinf(refimg)) = 0;

% bound image intensity with colormap bounds
img_cm = (img>cmax).*cmax + (img<cmin).*cmin + (img<cmax & img>cmin).*img;
%disp([ min(min(img_cm)) max(max(img_cm))])

% assign colormap index
img_cmidx = fix((img_cm-cmin)/(cmax-cmin)*cm_length)+1;

% transform to RGB
img_rgb = real(ind2rgb(img_cmidx,cm));

% transform to HSV
img_hsv = rgb2hsv(img_rgb);
% img_hsv(isnan(img_hsv))=0;  %is this needed?

img_brightness = img_hsv(:,:,3);
% img_brightness(isnan(Ratio)) = 0;  %one may want to set the brighness to 0 for ratios that are NaNs

%rescale reference image so that brightness from min_cutoff to cutoff ranges from 0 to 1
tmp = refimg;
tmp = (tmp>maxcutoff).*maxcutoff+(tmp<maxcutoff).*tmp; %ceiling
tmp = (tmp>maxcutoff).*maxcutoff+(tmp<mincutoff).*mincutoff+(tmp<=maxcutoff & tmp>=mincutoff).*tmp; %ceiling and floor
tmp = (tmp-mincutoff)./(maxcutoff-mincutoff); % rescale
% tmp(isnan(Ratio)) = 0; %one may want to set the brighness to 0 for ratios that are NaNs

% Apply brightness correction
img_hsv = cat(3,img_hsv(:,:,1:2),tmp);

% Transform back to RGB
img2 = hsv2rgb(img_hsv);


end


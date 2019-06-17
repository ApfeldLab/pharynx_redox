function imgStackRGB = imageStackToRGB(imgStack, map, range)
% IMAGESTACKTORGB Return an [height, width, color, frame] matrix, where the third dimension
% encodes RGB, and z gives image index
%
% ARGUMENTS:
%   imgStack
% a [height, width, frame] matrix
% 
%   map 
% a colormap, typically a [c x 3] matrix, where c can be up to 256
%
%   range 
% a two-element vector [min max], specifying the minimum and maximum pixel
% values that should be mapped to the bottom and top of the colormap 
% respectively. If none are supplied, this function sets it to be the 
% minimum and maximum value of the entire image stack.

    if nargin < 3
        range = [min(imgStack, [], [1 2 3]) max(imgStack, [], [1 2 3])];
    end
    
    sz = size(imgStack);
    
    imgStackInd = gray2ind(mat2gray(imgStack, range), size(map,1));
    imgStackRGB = zeros(sz(1), sz(2), 3, sz(3));

    for i=1:sz(3)
        imgStackRGB(:,:,:,i) = ind2rgb(imgStackInd(:,:,i), map);
    end
end
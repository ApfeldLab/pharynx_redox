function seg = segmentPharynx(imStack, useAC)
    seg = zeros(size(imStack));

    acSE = strel('square', 3);
    
    f = waitbar(0,'1','Name','Segmenting Images...',...
        'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
    setappdata(f,'canceling',0);

    steps = size(imStack, 3);

    h = fspecial('unsharp');
    for i=1:steps
        % Check for clicked Cancel button
        if getappdata(f,'canceling')
            break
        end

         % Update waitbar and message
        waitbar(i/steps,f, sprintf('Segmenting Animal %d', i));
%         mask = imfill(imdilate(imclose(edge(imStack(:,:,i), 'canny', [0.005 .1]), cannyse), cannyse), 'holes');
        mask = imfill(imdilate(imclose(edge(imStack(:,:,i), 'sobel'), strel('diamond',2)), strel('diamond',2)), 'holes');
        mask = bwpropfilt(logical(mask), 'Area', 1);
        
        if useAC
            unsharped = imfilter(imStack(:,:,i), h);
            seg(:,:,i) = imdilate(activecontour(unsharped, mask, 15, 'Chan-Vese', 'SmoothFactor', 3), acSE);
            seg(:,:,i) = bwpropfilt(logical(seg(:,:,i)), 'Area', 1);
        else
            seg(:,:,i) = mask;
        end
    end
    delete(f);

end
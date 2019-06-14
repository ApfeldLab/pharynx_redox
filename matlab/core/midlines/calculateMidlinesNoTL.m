function midlines = calculateMidlinesNoTL(fluorescenceSegStack)
    nAnimals = size(fluorescenceSegStack, 3);
    midlines = cell(nAnimals, 1);
    
    imHeight = size(fluorescenceSegStack, 1);
    imWidth  = size(fluorescenceSegStack, 2);
    
    textprogressbar('Calculating Midlines: ');
    for i=1:nAnimals
        textprogressbar(100 * (i/nAnimals));

        [y, x] = find(fluorescenceSegStack(:,:,i));
        
        % Add boundary points to constrain ends
        x = vertcat(x,0,imWidth);
        y = vertcat(y,imHeight/2,imHeight/2);
        
        [fitresult, ~] = fit(x, y, 'smoothingspline', 'SmoothingParam', 0.001);
        midlines(i) = {fitresult};
    end
    textprogressbar(' done');
end
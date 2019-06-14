function splitImages = splitImageStack(allImages, imagingStrategyString)
% TODO Documentation
    nImages = size(allImages, 3);
    lambdas = split(imagingStrategyString, "/");
    nLambdas = size(lambdas, 1);
    
    % Sanity Check
    assert(mod(nImages, nLambdas) == 0, ...
        sprintf("Number of frames in image stack (%d) is not divisible by number of wavelengths specified (%d).", nImages, nLambdas));

    splitImages = struct;

    lambdaCounter = containers.Map(cellstr(unique(lambdas)), num2cell(zeros(size(unique(lambdas), 1), 1)));

    for i = 1:nLambdas
        lambda = lambdas(i);
        lambdaCounter(lambda) = lambdaCounter(lambda) + 1;
        varname = strcat("im", lambda, "_", num2str(lambdaCounter(lambda)));
        splitImages.(varname) = allImages(:,:,i:size(lambdas,1):end);
    end
end
function data = analyze_image_stack(image_path, varargin)
% ANALYZE_IMAGE_STACK  Analyze pictures of the C. elegans pharynx for 
% analysis of reduction potentials 

persistent p
if isempty(p)
    p = inputParser;
    checknum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%     check_frame_specifier = @(x) strcmp(x, 'all') || (checknum(x) && isinteger(x));
    
    [default_output_dir, ~, ~] = fileparts(image_path);
    
    check_mask_threshold = @(x) strcmp(x, 'edge') || checknum(x);
    
    addRequired(p, 'ImagingScheme', @isstring);
    addOptional(p, 'OutputDir', default_output_dir, @(x) isstring(x) && isfolder(x));
    addOptional(p, 'SubtractMedians', false, @islogical);
    addOptional(p, 'RotateBeforeMeasurement', true, @islogical);
    addOptional(p, 'ProfileThreshold', 2000, checknum);
    addOptional(p, 'ChannelRegister', false, @islogical);
    addOptional(p, 'NPointsToSample', 1000, checknum);
    addOptional(p, 'MaskThreshold', 'edge', check_mask_threshold);
    addOptional(p, 'UseTLForMidlines', false, @islogical);
    addOptional(p, 'MeasurementInterpolationMethod', 'BILINEAR');
    addOptional(p, 'RotationInterpolationMethod', 'BILINEAR');
    addOptional(p, 'LengthNormalizeInterpMethod', 'linear');
end
parse(p, varargin{:});

data = struct;

data.images = splitImageStack(loadTiffImageStack(image_path), p.Results.ImagingScheme);
imageStackNames = fieldnames(data.images);
nonTLStackIdx = getNonTLidx(imageStackNames);

% Subtract Medians
if p.Results.SubtractMedians
    for i=size(imageStackNames, 1)
        data.images.(imageStackNames{i}) = subtractMedians(data.images.(imageStackNames{i}));
    end
end

% Segment
if strcmp(p.Results.MaskThreshold, 'edge')
    mask_thresh = 0; 
else
    mask_thresh = p.Results.MaskThreshold;
end

for i=1:size(imageStackNames, 1)
    stackName = imageStackNames{i};
    % TODO: segment TL differently than FL
    data.segImages.(stackName) = segmentPharynx(data.images.(stackName), 0, mask_thresh);
end

% Rotate Pharynxes
% TODO: rotate TL with FL
for i=1:length(nonTLStackIdx)
    idx = nonTLStackIdx(i);
    stackName = imageStackNames{idx};
    if p.Results.RotateBeforeMeasurement
        [rot_FL, rot_seg] = rotatePharynx(data.images.(stackName), data.segImages.(stackName));
        data.rot_FL.(stackName) = rot_FL;
        data.rot_seg.(stackName) = rot_seg;
    end
end

if p.Results.RotateBeforeMeasurement
    measurement_FL = data.rot_FL;
    measurement_seg = data.rot_seg;
else
    measurement_FL = data.images;
    measurement_seg = data.segImages;
end

% Calculate Midlines
framesForMidlines = getNonTLidx(imageStackNames);

for i=1:length(framesForMidlines)
    idx = framesForMidlines(i);
    stackName = imageStackNames{idx};
    if p.Results.UseTLForMidlines
        tlStack = data.images.(getCorrespondingTLName(stackName));
        
        data.midlines.(stackName) = calculateMidlines(tlStack, measurement_seg.(stackName), 0);
    else
        data.midlines.(stackName) = calculateMidlinesNoTL(measurement_seg.(stackName));
    end
end

% Measure Under Midlines
for i=1:length(framesForMidlines)
    idx = framesForMidlines(i);
    stackName = imageStackNames{idx};
    data.rawIntensities.(stackName) = measureAndTrim(...
        measurement_FL.(stackName), ...
        data.midlines.(stackName), ...
        p.Results.ProfileThreshold, ...
        p.Results.NPointsToSample, ...
        p.Results.MeasurementInterpolationMethod, ...
        p.Results.LengthNormalizeInterpMethod);
end

end


% Helper Functions
function stackName = getCorrespondingTLName(flName, stackNames)
    stackName = 'imTL_1';
end
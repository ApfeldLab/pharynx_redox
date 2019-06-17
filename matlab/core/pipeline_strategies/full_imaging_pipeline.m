analyze_image_stack(fullfile("data", "2017_02_22-HD233_SAY47.tif"),... 
    "TL/470/410/470/410");
"470/470/470/410/470" -> "470_1/470_2/470_3/410_1/470_4"


function analyze_image_stack(image_path, varargin)

persistent p
if isempty(p)
    p = inputParser;
    checknum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    check_frame_specifier = @(x) strcmp(x, 'all') || (checknum(x) && isinteger(x));
    
    [default_output_dir, ~, ~] = fileparts(image_path);
    
    check_mask_threshold = @(x) strcmp(x, 'edge') || checknum(x);
    
    addRequired(p, 'ImagingScheme', @isstring);
    addOptional(p, 'OutputDir', default_output_dir, @(x) isstring(x) && isfolder(x));
    addOptional(p, 'FramesForMidlines', 'all', check_frame_specifier); % [ 1 2 4 ]
    addOptional(p, 'FramesForMasks', 'all', check_frame_specifier);
    addOptional(p, 'SubtractMedians', false, @islogical);
    addOptional(p, 'RotateBeforeMeasurement', true, @islogical);
    addOptional(p, 'ProfileThreshold', 2000, checknum);
    addOptional(p, 'ChannelRegister', false, @islogical);
    addOptional(p, 'NPointsToSample', 1000, checknum);
    addOptional(p, 'MaskThreshold', 'edge', check_mask_threshold);
    addOptional(p, 'UseTLForMidlines', false, @islogical);
end
parse(p, varargin{:});

splitImageStacks = splitImageStack(loadTiffImageStack(image_path), p.Results.ImagingScheme);
imageStackNames = fieldnames(splitImageStacks);

% Subtract Medians
if p.Results.SubtractMedians
    for i=size(imageStackNames, 1)
        splitImageStacks.(imageStackNames{i}) = subtractMedians(splitImageStacks.(imageStackNames{i}));
    end
end

% Segment
if strcmp(p.Results.MaskThreshold, 'edge')
    mask_thresh = 0; 
else
    mask_thresh = p.Results.MaskThreshold;
end

segImages = struct;

if strcmp(p.Results.FramesForMasks, 'all')
    framesToSegment = 1:size(imageStackNames, 1);
else
    framesToSegment = p.Results.FramesForMasks;
end

for i=1:size(framesToSegment, 2)
    idx = framesToSegment(i);
    stackName = imageStackNames{idx};
    % TODO: segment TL differently than FL
    segImages.(stackName) = segmentPharynx(splitImageStacks.(stackName), 0, mask_thresh);
end

% Calculate Midlines
midlines = struct;

if strcmp(p.Results.FramesForMasks, 'all')
    framesForMidlines = 1:size(imageStackNames, 1);
else
    framesForMidlines = p.Results.FramesForMidlines;
end

for i=1:size(framesForMidlines, 2)
    idx = framesForMidlines(i);
    stackName = imageStackNames{idx};
    if p.Results.UseTLForMidlines
        tlStack = splitImageStacks.(getCorrespondingTLName(stackName));
        
        % TODO: this is broken if only a subset of stacks have been segmented
        midlines.(stackName) = calculateMidlines(tlStack, segImages.(stackName), 0);
    else
        midlines.(stackName) = calculateMidlinesNoTL(splitImageStacks.(stackName));
    end
end

% Measure Under Midlines
for i=1:size(imageStackNames,1)
    % If ~TL
    %   Get corresponding midline
    %       if there is a midline with same name as stackname, 
    %           that's it
    %       else
    %           find all midlines with same number
    
    % skip
end


end

% Helper Functions
function stackName = getCorrespondingTLName(flName, stackNames)
    
end
classdef Constants
    %CONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        regions = struct(...
            'pm3', [70 90], 'pm4', [59 65], ...
            'pm5', [30 50], 'pm6', [15 20], ...
            'pm7', [5 10],  'medial_axis', [3 97]);
        
        lambda = 0.0891;
    end
end
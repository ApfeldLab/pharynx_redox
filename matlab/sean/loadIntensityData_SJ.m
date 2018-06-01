function sqData = loadIntensityData_SJ(filepath, t)
%loadIntensityData_SJ Load and clean intensity data from MultiMeasure
%output of ImageJ
%   ARGUMENTS:
%   t = threshold
    
    % Ignore Headers
    data = csvread(filepath, 1);
    
    % Remove X Values
    yData = data;
%     yData = data(:, 2:2:end);
    
    % Threshold
%     for i=1:size(yData,2)
%         yData(1:find(yData(:,i)>t, 1, 'first'), i) = 0;
%         yData(find(yData(:,i)>t, 1, 'last'):end, i) = 0;
%     end
    
    % Length-Normalize
    sqData = ssquare(yData);
end


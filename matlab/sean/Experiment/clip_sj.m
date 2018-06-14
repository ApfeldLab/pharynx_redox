function clipped = clip_sj(data, t)
%CLIP_SJ Threshold, but not in the middle of the data
    clipped = data;
    for i=1:size(data,2)
        clipped(1:find(data(:,i)>t, 1, 'first'), i) = 0;
        clipped(find(data(:,i)>t, 1, 'last'):end, i) = 0;
    end
end
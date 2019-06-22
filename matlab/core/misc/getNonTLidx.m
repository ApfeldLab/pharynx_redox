function idx = getNonTLidx(stackNames)
    idx = {};
    for i=1:length(stackNames)
        name = stackNames{i};
        if ~contains(name, 'TL', 'IgnoreCase', true)
            idx{end + 1} = i; %#ok<AGROW>
        end
    end
    idx = cell2mat(idx);
    return
end
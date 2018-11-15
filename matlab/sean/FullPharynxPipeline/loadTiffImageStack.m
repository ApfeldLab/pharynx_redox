function imgStack = loadTiffImageStack(filepath)
    imTable = struct2table(tiffread2(filepath));
    imgStack = double(cat(3, imTable.data{:})); % dimensions: WxHxi; e.g. first image is ims(:, :, 1)
end
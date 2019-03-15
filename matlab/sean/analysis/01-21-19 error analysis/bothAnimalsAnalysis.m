rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);

allImgs = loadTiffImageStack(rawImgFilePath);

% Load Movement
m_sep = readtable("sean/analysis/01-21-19 error analysis/data/movement_separated.csv");

I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');


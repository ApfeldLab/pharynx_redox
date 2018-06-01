roiManager("reset");
dirName = getInfo("image.directory");
runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
segTitle = getTitle;
setThreshold (25, 60000);

run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=�m");// 1px=2.58�m for 4x4 ; 1px=1.29�m for 2x2
run("Set Measurements...", "area centroid center perimeter bounding fit shape feret's median skewness kurtosis scientific redirect=None decimal=7");
run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display add stack");
run("Clear Results");
roiManager("deselect");
roiManager("measure");

selectWindow("Results");
save(dirName + "morph.txt");
run("Close");
selectWindow(segTitle);
close();
roiManager("reset");
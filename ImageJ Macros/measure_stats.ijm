dirName = getInfo("image.directory");
imTitles = getList("image.titles");

runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
segTitle = getTitle;
setThreshold (25, 60000);
run("Set Measurements...", "area centroid center perimeter bounding fit shape feret's median skewness kurtosis scientific redirect=None decimal=7");
run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display add stack");
run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7"); 
for (j = 0; j < nSlices; j++) {
	setSlice(j+1);
	run("Measure");
}
selectWindow("Results");
save(dirName + "measurements.txt");
run("Close");
selectWindow(segTitle);
close();
roiManager("reset");
dirName = getInfo("image.directory");

run("Clear Results");
// Put the coordinates of all ROI lines in a Results window
title = getTitle;
for (i = 0; i < roiManager("count"); i++) {
	roiManager("select", i);
	getSelectionCoordinates(x, y);
	for (j = 0; j < lengthOf(x); j++) {
		setResult("X"+i, j, x[j]);
		setResult("Y"+i, j, y[j]); 	
	}
}
updateResults;

baseName = split(getTitle, ".");
baseName = baseName[0];
savePath =  dirName + baseName + "-polyCoords.csv";
saveAs("results", savePath);
selectWindow("Results");
run("Close");
showMessage("Saved Polyline Coordinates to: " + savePath);
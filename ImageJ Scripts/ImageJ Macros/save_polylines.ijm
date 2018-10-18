run("Set Scale...", "distance=0 global");
// With the image stack selected, make sure a polyline has already been drawn, then run this macro.
// It will save the polyline coordinates and measurements with the correct formats & names.
dirName = getInfo("image.directory");
title = getTitle;
baseTitle = split(title, '.');
baseTitle = baseTitle[0];

// SAVE COORDS
for (i = 0; i < roiManager("count"); i++) {
	roiManager("select", i);
	getSelectionCoordinates(x, y);
	for (j = 0; j < lengthOf(x); j++) {
		setResult("X"+i, j, x[j]);
		setResult("Y"+i, j, y[j]);
	}
}
updateResults;
selectWindow("Results");
save(dirName + baseTitle + "_coords.txt");
run("Close");

// SAVE MEASUREMENTS
run("Clear Results");
for (i = 0; i < roiManager("count"); i++) {
	roiManager("select", i);
	profile = getProfile();
	for (j=0; j<profile.length; j++) {
		setResult(""+i, j, profile[j]);
	}
}

updateResults;
selectWindow("Results");
save(dirName + baseTitle + "_intensities.txt");
run("Close");
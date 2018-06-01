dirName = getInfo("image.directory");
title = getTitle;
baseTitle = split(title, '.');
baseTitle = baseTitle[0];

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
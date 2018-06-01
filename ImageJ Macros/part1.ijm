run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=cm global");

dirName = getInfo("image.directory");

runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/Split_Stack.ijm", "stack");
titles = getList("image.titles");
Array.print(titles);
// Initial Measurements
run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7");

for (i = 0; i < lengthOf(titles);  i++) {
	selectWindow(titles[i]);
	baseTitle = split(titles[i], '.');
	baseTitle = baseTitle[0];
	
	for (j = 0; j < nSlices; j++) {
		setSlice(j+1);
		run("Measure");
	}
	selectWindow("Results");
	save(dirName + baseTitle + ".txt");
	run("Close");
}

// Subtract Mode
run("Set Measurements...", "modal redirect=None decimal=7");
for (i = 0; i < lengthOf(titles); i++) {
	selectWindow(titles[i]);
	baseTitle = split(titles[i], '.');
	baseTitle = baseTitle[0];
	
	run("Duplicate...", "duplicate");
	for (j = 0; j < nSlices; j++) {
		setSlice(j+1);
		run("Measure");
		Background = getResult("Mode");
        run("Subtract...", "value="+Background);
	}
	rename(baseTitle + "_subMode.tif");
	save(dirName + baseTitle + "_subMode.tif");
}

selectWindow("Results");
run("Close");

// Measure segmented pharynxes of Mode-Subtracted Images
run("Set Measurements...", "area mean standard modal min centroid center bounding fit redirect=None decimal=7");
for (i = 0; i < lengthOf(titles); i++) {
	baseTitle = split(titles[i], '.');
	baseTitle = baseTitle[0];
	if (!endsWith(baseTitle, "TL")) {
		roiManager("reset");
		run("Clear Results");
		selectWindow(baseTitle + "_subMode.tif");
		runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
		segTitle = getTitle;
		setThreshold (25, 60000);
		run("Analyze Particles...", "size=200-Infinity circularity=0.00-1.00 show=Nothing display add stack");
		roiManager("Measure");
		selectWindow("Results");
		save(dirName + baseTitle + "pharynxMeasurements.txt");
		run("Close");

		selectWindow(segTitle);
		close();
	}	
}

showMessage("Substacks Created, Measurements Saved");
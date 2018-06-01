dirName = getInfo("image.directory");
imTitles = getList("image.titles");


Dialog.create("Measure Modes");
Dialog.addMessage("Please which image stacks to measure the mode.")
defaults = newArray(lengthOf(imTitles));
Array.fill(defaults, true);
Dialog.addCheckboxGroup(lengthOf(imTitles), 1, imTitles, defaults);
Dialog.show();

run("Set Measurements...", "modal redirect=None decimal=7");

shouldProcessIdx = newArray(lengthOf(imTitles));
for (i = 0; i < lengthOf(imTitles); i++) {
	shouldProcessIdx[i] = Dialog.getCheckbox();
}

for (i = 0; i < lengthOf(shouldProcessIdx); i++) {
	if (shouldProcessIdx[i] == 1) {
		selectWindow(imTitles[i]);
		run("Clear Results");
		for (j=1; j < nSlices - 1; j++) {
			setSlice(j);
			run("Measure");
		}
		selectWindow("Results");
		baseName = split(imTitles[i], '.');
		save(dirName + baseName[0] + "-modes.txt");
	}
}
selectWindow("Results");
run("Close");
showMessage("Saved Modes");
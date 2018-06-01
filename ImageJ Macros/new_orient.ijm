//////////////////////////////////////////////////////////////
// Setup
//////////////////////////////////////////////////////////////
dirName = getInfo("image.directory");
imTitles = getList("image.titles");

Dialog.create("Rotate Images");
Dialog.addMessage("Please select the image stack you would like the rotation to be based on.")
Dialog.addChoice("Rotation Image Stack", imTitles);


Dialog.addMessage("Please select the image stacks you would like to rotate/crop.")
defaults = newArray(lengthOf(imTitles));
Array.fill(defaults, true);
Dialog.addCheckboxGroup(lengthOf(imTitles), 1, imTitles, defaults);
Dialog.show();

rotStackTitle = Dialog.getChoice();
shouldProcessIdx = newArray(lengthOf(imTitles));
for (i = 0; i < lengthOf(imTitles); i++) {
	shouldProcessIdx[i] = Dialog.getCheckbox();
}

//////////////////////////////////////////////////////////////
// Calculate Angles and Centers of Rotation
//////////////////////////////////////////////////////////////
selectWindow(rotStackTitle);

// Segment images
runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
segTitle = getTitle;
setThreshold (25, 60000);

run("Clear Results");
run("Set Measurements...", "fit centroid shape");
run("Analyze Particles...", "size=200-Infinity pixel circularity=0.00-1 show=Nothing display add stack");

roiManager("show none");
roiManager("reset");
selectWindow(segTitle);
close();

//////////////////////////////////////////////////////////////
// Do the Rotation and Cropping for each stack
//////////////////////////////////////////////////////////////
for (i = 0; i < lengthOf(shouldProcessIdx); i++) {
	if (shouldProcessIdx[i] == 1) {
		rotateStack(imTitles[i]);
		cropStack("PA-" + imTitles[i]);
	}
}

//////////////////////////////////////////////////////////////
// Posterior/Anterior Alignment for each stack
//////////////////////////////////////////////////////////////
run("Clear Results");
run("Set Measurements...", "area shape redirect=None decimal=7");

selectWindow("PA-"+rotStackTitle);
runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
segTitle = getTitle;
setThreshold(256, 60000);

// Build up an array that tells us whether or not to flip each image in the stack
shouldFlip = newArray(nSlices);
for (i = 0; i < nSlices; i++) {	
	run("Clear Results");
	setSlice(i + 1);
    makeRectangle(5,9,25,25);
    run("Analyze Particles...", "size=10-Infinity pixel circularity=0.00-1.00 show=Nothing display");	run("Select None");	

    makeRectangle(42,9,25,25);
    run("Analyze Particles...", "size=10-Infinity pixel circularity=0.00-1.00 show=Nothing display");	run("Select None");	

    leftEnd = getResult("Round", 0);
    rightEnd = getResult("Round", 1);	
    if (leftEnd < rightEnd){
        shouldFlip[i] = true;
   	} else {
   		shouldFlip[i] = false;
   	}
}
resetThreshold();

// Do the actual flipping
for (i = 0; i < lengthOf(shouldProcessIdx); i++) {
	if (shouldProcessIdx[i] == 1) {
		flipIfNecessary("PA-" + imTitles[i], shouldFlip);
	}
}

close(segTitle);

//////////////////////////////////////////////////////////////
// Save PA Images
//////////////////////////////////////////////////////////////
for (i = 0; i < lengthOf(shouldProcessIdx); i++) {
	if (shouldProcessIdx[i] == 1) {
		selectWindow("PA-" + imTitles[i]);
		save(dirName + "PA-" + imTitles[i]);
	}
}

showMessage("Saved PA Alignments\nMake sure to go through the images, check alignments, and resave if necessary.");


//////////////////////////////////////////////////////////////
// Helper functions
//////////////////////////////////////////////////////////////
function rotateStack(title) {
	selectWindow(title);
	width = getWidth;
	height = getHeight;
	
	xCenter = width / 2;
	yCenter = height / 2;
	run("Duplicate...", "duplicate title=PA-"+title);
	for (i = 0; i < nSlices; i++) {
		setSlice(i+1);
		
		//roiManager("Select", i);
		xOffset = xCenter - getResult("X", i);
		yOffset = yCenter - getResult("Y", i);
		Angle =  getResult("Angle", i);	
	
		run("Select None");
		run("Translate...", "x="+xOffset+" y="+yOffset+" interpolation=None slice");
		run("Rotate... ", "angle=Angle grid=1 interpolation=Bilinear slice");
	}
}

function cropStack(title) {
	CROP_H = 40;
	CROP_W = 75;
	
	crop_x = getWidth / 2 - CROP_W / 2;
	crop_y = getHeight / 2 - CROP_H / 2;

	selectWindow(title);
	makeRectangle(crop_x, crop_y, CROP_W, CROP_H);
	run("Crop");
}

function flipIfNecessary(title, shouldFlipArray) {
	selectWindow(title);
	for (j = 0; j < nSlices; j++) {
		setSlice(j + 1);
		if (shouldFlipArray[j]) {
			run("Flip Horizontally", "slice");
		}
	}
}

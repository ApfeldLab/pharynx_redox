// Given a stack of images, rotate and align them to the center

// Overall method:
// I => SegI
// Rotate(I) according to angle(segI) and CoM(SegI)
// Horizontally flip images so they are all Posterior -> Anterior

run("Set Scale...", "distance=0 global");
origTitle = getTitle;

// Segment images
runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
segTitle = getTitle;

// Get Angles and Centers of Rotation for the segmented images
setThreshold (25, 60000);
run("Clear Results");
run("Set Measurements...", "fit centroid shape");
run("Analyze Particles...", "size=200-Infinity pixel circularity=0.00-1 show=Nothing display add stack");
roiManager("show none");
roiManager("reset");
resetThreshold();
// Rotate Images
rotateStack(origTitle);
rotateStack(segTitle);

selectWindow("Results"); 
run("Close");

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

// Crop
CROP_H = 40;
CROP_W = 75;
CROP_X = getWidth / 2 - CROP_W / 2;
CROP_Y = getHeight / 2 - CROP_H / 2;
selectWindow("PA-"+origTitle);
makeRectangle(CROP_X, CROP_Y, CROP_W, CROP_H);
run("Crop");
selectWindow("PA-"+segTitle);
makeRectangle(CROP_X, CROP_Y, CROP_W, CROP_H);
run("Crop");

// Posterior -> Anterior Alignment
run("Clear Results");
run("Set Measurements...", "area shape redirect=None decimal=7");

selectWindow("PA-"+segTitle);
setThreshold(256, 60000);

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

flipIfNecessary("PA-"+origTitle, shouldFlip);
flipIfNecessary("PA-"+segTitle, shouldFlip);

function flipIfNecessary(title, shouldFlipArray) {
	selectWindow(title);
	for (j = 0; j < nSlices; j++) {
		setSlice(j + 1);
		if (shouldFlipArray[j]) {
			run("Flip Horizontally", "slice");
		}
	}
}
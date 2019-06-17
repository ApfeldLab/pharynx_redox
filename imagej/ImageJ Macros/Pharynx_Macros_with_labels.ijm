TITLES = newArray();



macro "Split [S]" { 
	run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=mm global");
	dirName = getInfo("image.directory"); 
	baseTitle = split(getTitle, "."); 
	baseTitle = baseTitle[0]; 
 
	Dialog.create("Split Channels"); 
	Dialog.addMessage("Enter the channels, and comma-separate if >1 (e.g. 2,3,5)") 
	Dialog.addString("TL Channel(s)", "1"); 
	Dialog.addString("470 Channel(s)", "2"); 
	Dialog.addString("410 Channel(s)", "3"); 
	Dialog.show(); 
 
	strListTL = Dialog.getString(); 
	strList470 = Dialog.getString(); 
	strList410 = Dialog.getString(); 
 
	listTL = split(strListTL, ","); 
	list470 = split(strList470, ","); 
	list410 = split(strList410, ","); 
 
	n_channels = lengthOf(listTL) + lengthOf(list470) + lengthOf(list410); 
 
	run("Deinterleave", "how=n_channels"); 
 
	for (i=0; i < lengthOf(listTL); i++) { 
		selectWindow(baseTitle + ".tif" + " #" + listTL[i]); 
		rename(baseTitle + "_TL.tif"); 
		save(dirName + getTitle); 
	} 
 
	for (i=0; i < lengthOf(list470); i++) { 
		selectWindow(baseTitle + ".tif" + " #" + list470[i]); 
		rename(baseTitle + "_470.tif"); 
		save(dirName + getTitle); 
	} 
 
	for (i=0; i < lengthOf(list410); i++) { 
		selectWindow(baseTitle + ".tif" + " #" + list410[i]); 
		rename(baseTitle + "_410.tif"); 
		save(dirName + getTitle); 
	} 
} 
 
macro "Make Masks [1]" { 
	run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=mm global");
	dirName = getInfo("image.directory"); 
	titles = getList("image.titles"); 
 
	for (i = 0; i < lengthOf(titles); i++) { 
		bName = baseName(titles[i]); 
		if (!endsWith(bName, "TL")) { 
			selectWindow(titles[i]); 
			makeMask(); 
			save(dirName + maskName(titles[i])); 
			close(); 
		} 
	} 
	showMessage("Saved Masks"); 
} 
 
macro "Initial Measurements and Mode Subtraction [2]" { 
	run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=mm global");
	dirName = getInfo("image.directory"); 
	titles = getList("image.titles"); 
	 
	// Initial Measurements 
	run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7"); 
	for (i = 0; i < lengthOf(titles);  i++) { 
		roiManager("reset");
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
	// Subtract Median 
	run("Set Measurements...", "median redirect=None decimal=7"); 
	for (i = 0; i < lengthOf(titles); i++) { 
		selectWindow(titles[i]); 
		baseTitle = split(titles[i], '.'); 
		baseTitle = baseTitle[0]; 
 
		run("Duplicate...", "duplicate"); 
		for (j = 0; j < nSlices; j++) { 
			setSlice(j+1); 
			run("Measure"); 
			Background = getResult("Median"); 
	        run("Subtract...", "value="+Background); 
		} 
		rename(baseTitle + "_subMed.tif"); 
		save(dirName + baseTitle + "_subMed.tif"); 
	} 
 
	selectWindow("Results"); 
	run("Close"); 
	// Measure segmented pharynxes of Mode-Subtracted Images 
 
	// Segment, and save segmented images 
	for (i = 0; i < lengthOf(titles); i++) { 
		baseTitle = split(titles[i], '.'); 
		baseTitle = baseTitle[0]; 
		if (!endsWith(baseTitle, "TL")) { 
			roiManager("reset"); 
			run("Clear Results"); 
			selectWindow(baseTitle + "_subMed.tif"); 
			threshold(); 
			save(dirName + segName(titles[i])); 
		} 
	} 
	 
	n_animals = nSlices;
	run("Set Measurements...", "area mean standard modal min centroid center bounding fit redirect=None decimal=7"); 
	for (i = 0; i < lengthOf(titles); i++) { 
		baseTitle = split(titles[i], '.'); 
		baseTitle = baseTitle[0]; 
		if (!endsWith(baseTitle, "TL")) { 
			selectWindow(baseTitle + "_subMed.tif"); 
			threshold(); 
			segTitle = getTitle; 
			setThreshold (250, 60000); 
			
			run("Clear Results");
			roiManager("reset");
			run("Analyze Particles...", "size=200-Infinity circularity=0.00-1.00 show=Nothing display add stack"); 
			selectWindow("Results"); 
			save(dirName + baseTitle + "_pMeasure.txt"); 
			
			run("Close"); 
			selectWindow(segTitle); 
			close(); 
			roiManager("reset"); 
		} 
	} 
 
	showMessage("Substacks Created, Measurements Saved"); 
} 
 
macro "Orient and Align [3]" { 
	run("Set Scale...", "distance=0 global");
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
	threshold();
	segTitle = getTitle; 
	setThreshold (250, 60000); 
 	roiManager("reset");
	run("Clear Results"); 
	run("Set Measurements...", "fit centroid shape"); 
	run("Analyze Particles...", "size=200-Infinity circularity=0.00-1 show=Nothing display add stack"); 
 
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
	threshold();
	segTitle = getTitle; 
	setThreshold(256, 60000); 
 
	// Build up an array that tells us whether or not to flip each image in the stack 
	shouldFlip = newArray(nSlices); 
	for (i = 0; i < nSlices; i++) { 
		run("Clear Results"); 
		setSlice(i + 1); 
	    //makeRectangle(5,9,25,25); 
	    makeRectangle(17,4,22,30);
	    run("Analyze Particles...", "size=10-Infinity pixel circularity=0.00-1.00 show=Nothing display");	run("Select None"); 
 
	    //makeRectangle(42,9,25,25); 
	    makeRectangle(50,4,26,29);
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
		CROP_W = 90; 
 
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
} 
 
macro "Measure Morphology [4]" { 
	run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=mm global");
	// Only need to run this on either 410 or 470 
	dirName = getInfo("image.directory"); 
	makeMask(); 
	run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=�m");// 1px=2.58 m for 4x4 ; 1px=1.29 m for 2x2 
	run("Set Measurements...", "area centroid center perimeter bounding fit shape feret's median skewness kurtosis scientific redirect=None decimal=7"); 
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display add stack"); 
	close(); 
	selectWindow("Results"); 
	save(dirName + "morph.txt"); 
	run("Close"); 
	roiManager("reset");
} 
 
macro "Draw Polyline [P]" { 
	run("Set Scale...", "distance=0 global");
	// Select Image window to draw on, then run macro 
	width = getWidth; 
	height = getHeight; 
	title = getTitle(); 

	run("Line Width...", "line=2"); 
	 

	Dialog.create("Select Threshold PA");
	titles = getList("image.titles");
	Dialog.addChoice("Threshold PA", titles);
	Dialog.show();
	maskPaTitle = Dialog.getChoice();

	selectWindow(maskPaTitle);
	
	for (j = 0; j < nSlices; j++) { 
		setSlice(j+1); 
	 
		polyline_slice(); 
	}
	
	//selectWindow("Results"); 
	//run("Close");
	 
	function polyline_slice() { 
		width = getWidth; 
		height = getHeight; 
		title = getTitle(); 
		run("Clear Results"); 
		 
		// Get Center of Masses 
		run("Set Measurements...", "center mean"); 
		for (i = 0; i < width; i++) { 
			makeRectangle(i, 0, 1, height); 
			run("Measure"); 
		} 
	 
		leftBound = 0; 
		for (i = 0; i < nResults; i++) { 
			if (getResult("Mean", i) != 0.0) { 
				leftBound = i; 
				break; 
			} 
		} 
	 
		rightBound = 0; 
		for (i = nResults - 1; i >= 0; i--) { 
			if (getResult("Mean", i) != 0.0) { 
				rightBound = i; 
				break; 
			} 
		} 
		 
		Xs = newArray(nResults); 
		Ys = newArray(nResults); 
		 
		for (i = 0; i < nResults; i++) { 
			Xs[i] = i; 
			Ys[i] = getResult("YM", i); 
		} 
		 
		XTrimmed = Array.slice(Xs, leftBound + 7, rightBound); 
		YTrimmed = Array.slice(Ys, leftBound + 7, rightBound); 
	 
		pad = 6; 
		lineWidth = rightBound - leftBound + (pad * 2);
		Fit.doFit("3rd Degree Polynomial", XTrimmed, YTrimmed); 
	 
		XFit = newArray(width); 
		YFit = newArray(width); 
	 
		nPts = 6; 
		XFit = newArray(nPts+1); 
		YFit = newArray(nPts+1); 
		 
		for (i = 0; i <= nPts; i++) { 
			XFit[i] = (leftBound - pad) + ((lineWidth / nPts) * i); 
			YFit[i] = minOf(Fit.f(XFit[i]), height - 10); 
		} 
	 
		makeSelection(6, XFit, YFit); 
		run("Fit Spline"); 
		roiManager("add"); 
		reset; 
	} 
}
 
macro "Save Polyline [Q]" { 
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
} 
 
macro "Reset Environment [r]" { 
	// Close all images 
	while (nImages>0) {  
    	selectImage(nImages);  
        close();  
    } 
    run("Clear Results"); 
    selectWindow("Results"); 
    run("Close"); 
    roiManager("deselect"); 
    roiManager("delete"); 
    selectWindow("ROI Manager"); 
    run("Close"); 
    run("Main Window [return]"); 
} 
 
macro "Test Threshold [t]" { 
	// Runs the threshold on the image stack. 
	threshold();
} 

macro "merge" {
	run("Merge Channels...", "c1=[Mask of vab1_2do_05_25_mvmt_470-2.tif] c4=vab1_2do_05_25_mvmt_TL.tif create keep");
}

// HELPER FUNCTIONS 

// Create a new window, with a binary mask [0, 255]
// called "Mask of [original title]"
function threshold() { 
	origTitle = getTitle;
	run("Duplicate...", "duplicate");
	dupID = getImageID;
	AUTOFLORESCENCE_THRESH = 1000;
	PHARYNX_THRESH = 1000;
	run("Find Edges", "stack");
	setThreshold(1000,100000);
	run("Analyze Particles...", "size=50-Infinity show=Masks exclude stack");
	run("Open", "stack");
	selectImage(dupID);
	close();
}
 
function makeMask() { 
	// Assumes that the window to run this on is open 
	threshold();
	//segTitle = getTitle; 
	//setThreshold(256, 600000); 
	//run("Analyze Particles...", "size=200-Infinity pixel show=Masks stack"); 
	//maskTitle = getTitle; 
	//close(segTitle); 
}
 
// FILENAME FUNCTIONS 
function baseName(windowName) { 
	// Removes the file extension from the window name 
	// foobar.tif -> foobar 
	// foobar.txt -> foobar 
	baseTitle = split(windowName, '.'); 
	bName = baseTitle[0]; 
	return bName; 
} 
 
function segName(windowName) { 
	bName = baseName(windowName); 
	return bName + "-SEG.tif"; 
} 
 
function maskName(windowName) { 
	bName = baseName(windowName); 
	return bName + "-MASK.tif"; 
} 
 
function alignedName(windowName) { 
	bName = baseName(windowName); 
	return "PA-" + bName + ".tif"; 
} 
 
function subModeName(windowName) { 
	bName = baseName(windowName); 
	return bName + "_nomode.tif"; 
} 
 
function pharynxMeasurementsName(windowName) { 
	bName = baseName(windowName); 
	return bName + "_pMeasure.txt"; 
} 
 
function totalImageMeasurementsName(windowName) { 
	bName = baseName(windowName); 
	return bName + ".txt"; 
} 
 
function polyCoordsName(windowName) { 
	bName = baseName(windowName); 
	return bName + "_coords.txt"; 
} 
 
function polyMeasureName(windowName) { 
	bName = baseName(windowName); 
	return bName + "_intensities.txt"; 
} 
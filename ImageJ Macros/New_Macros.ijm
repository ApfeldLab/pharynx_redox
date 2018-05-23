//open("/Users/sean/code/wormAnalysis/imageJ_Macro_Playground/HD233_High_Movement_03_16_by5.tif");

N_CHANNELS = 3;

macro "RotateAlign_PA" {
	// Split Stack into 5 Channels
    frames = nSlices/N_CHANNELS;
	run("Stack to Hyperstack...", "order=xyczt(default) channels=N_CHANNELS slices=1 frames=frames display=Color");
    run("Split Channels");
    titles = getImageWindowTitles();
    
    //subtractModes();
    
    //RoiForStack(titles[N_CHANNELS - 1]);
    RoiForStack(titles[1]);

    for (i = 0; i < titles.length; i++) {
    	rotateStack(titles[i]);
    }
    arrangeTranslatedWindows();
    alignPA();
    //divideTranslated();
}

function subtractModes() {
	run("Set Measurements...", "modal redirect=None decimal=7");
	for (i = 0; i < titles.length; i++) {
		selectWindow(titles[i]);
		for (j = 0; j < nSlices; j++) {
			setSlice(j+1);
			run("Clear Results");
			run("Measure");
	        Background = getResult("Mode");
	        run("Subtract...", "value=Background");
	        resetMinAndMax();
		}
	}
}

function getImageWindowTitles() {
	return getList("image.titles");
}

function getImageIDs() {
	IDs = newArray(nImages()); 
  	for (i=1; i<=nImages(); i++) { 
    	selectImage(i); 
    	IDs[i-1] = getImageID();
    }
    return IDs;
}

function RoiForStack(title) {
	selectWindow(title);
	setSlice(1);
	
	//setThreshold (1000, 60000);
	setAutoThreshold("Li dark"); 
	run("Clear Results");
	run("Set Measurements...", "fit centroid shape");
	run("Analyze Particles...", "size=200-Infinity circularity=0.00-1 show=Nothing display add stack");
	roiManager("show none");
	resetThreshold();

	// Remove ROIs that are too "blobby", which results in an ROI around the body not pharynx
	// Have to iterate backwards so as not to mess up the order as we delete
	for (i = 0; i < nResults; i++) {
		roundness = getResult("Round", i);
		if (roundness > 0.4) {
			roiManager("select", i);
		}
	}
	roiManager("delete");
	setSlice(1);
}

function rotateStack(title) {
	selectWindow(title);
	width = getWidth;
	height = getHeight;

	xCenter = width / 2;
	yCenter = height / 2;
	
	selectWindow(title);

	//translateName = getTranslatedName(title);
	
	//newImage(translateName, bitDepth, width, height, nResults);

	//run("Duplicate...", "duplicate title="+translateName);
	run("Duplicate...", "duplicate title="+translateName);

    for (i = 0; i < nResults; i++) {
    	setSlice(i+1);

		roiManager("Select", i);
		xOffset = xCenter - getResult("X", i);
		yOffset = yCenter - getResult("Y", i);
		Angle =  getResult("Angle", i);	

		run("Select None");
		run("Translate...", "x="+xOffset+" y="+yOffset+" interpolation=None slice");
		run("Rotate... ", "angle=Angle grid=1 interpolation=Bilinear slice");
    }

    R_WIDTH = 70;
    R_HEIGHT = 35;
    makeRectangle((width / 2) - floor(R_WIDTH / 2), (height / 2) - floor(R_HEIGHT / 2), R_WIDTH, R_HEIGHT);
	run("Crop");
	run("Select None");
	resetMinAndMax();
}

function arrangeTranslatedWindows() {
	Xs = newArray(278, 520, 278, 520, 278);
	Ys = newArray(273, 273, 450, 450, 627);
	
	for (i = 0; i < titles.length; i++) {
		if (isOpen("TRANSLATED-" + titles[i])); {
			selectWindow("TRANSLATED-" + titles[i]);
			setLocation(Xs[i], Ys[i], 300, 176);
		}
	}
}

function alignPA() {
	selectWindow(getTranslatedName(titles[N_CHANNELS-2]));
	run("Set Measurements...", "area shape redirect=None decimal=7");
	//setThreshold(1200, 60000);
	setAutoThreshold("Li dark"); 
	shouldFlip = newArray(nSlices);
	for (i = 0; i < nSlices; i++) {	
		run("Clear Results");
		setSlice(i + 1);
	    //makeRectangle(0, 29, 37, 45);
	    makeRectangle(6, 6, 17, 22); // for old
	    run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing display");	run("Select None");	
	
	    //makeRectangle(62, 29, 37, 45);
	    makeRectangle(46, 6, 17, 22); // for old
	    run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing display");	run("Select None");	
	
	    leftEnd = getResult("Round", 0);
	    rightEnd = getResult("Round", 1);	
	    if (leftEnd < rightEnd){
	        shouldFlip[i] = true;
	   	} else {
	   		shouldFlip[i] = false;
	   	}
	}
	resetThreshold();

	for (i = 0; i < titles.length; i++) {
		selectWindow(getTranslatedName(titles[i]));
		for (j = 0; j < nSlices; j++) {
			setSlice(j + 1);
			if (shouldFlip[j]) {
				run("Flip Horizontally", "slice");
			}
		}
	}
}

function divideTranslated() {
	// 410_1 / 470_1
	//divideByTitle(getTranslatedName(titles[1]), getTranslatedName(titles[2]));
	//rename("R410_470_1_"+substring(titles[1], 3));
	// 410_2 / 470_2
	//divideByTitle(getTranslatedName(titles[3]), getTranslatedName(titles[4]));
	//rename("R410_470_2_"+substring(titles[1], 3));

	// 410_1 / 410_2
	divideByTitle(getTranslatedName(titles[1]), getTranslatedName(titles[3]));
	rename("R410_410_"+substring(titles[1], 3));
	// 470_1 / 470_2
	divideByTitle(getTranslatedName(titles[2]), getTranslatedName(titles[4]));
	rename("R470_470_"+substring(titles[1], 3));
}

function divideByTitle(title1, title2) {
	selectWindow(title1);
	id1 = getImageID();
	selectWindow(title2);
	id2 = getImageID();
	
	imageCalculator("divide create 32-bit stack", id1, id2);
}

function getTranslatedName(title) {
	return "TRANSLATED-" + title;
}

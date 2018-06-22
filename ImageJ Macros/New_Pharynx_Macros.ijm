resetMinAndMax();

DIR = getInfo("image.directory");
TITLES = newArray();

macro "Split [0]" {
	dirName = getInfo("image.directory");
	baseTitle = split(getTitle, ".");
	baseTitle = baseTitle[0];

	Dialog.create("Split Channels");
	Dialog.addMessage("Enter the channels index")
	Dialog.addString("TL Channel", "1");
	Dialog.addString("470 Channel", "2");
	Dialog.addString("410 Channel", "3");
	Dialog.show();

	channel_tl = Dialog.getString();
	channel_470 = Dialog.getString();
	channel_410 = Dialog.getString();

	n_channels = 3;

	run("Deinterleave", "how=n_channels");

	// Rename and Save Stacks
	selectWindow(baseTitle + ".tif" + " #" + channel_tl);
	rename(baseTitle + "_TL.tif");
	resetMinAndMax();
	save(dirName + getTitle);

	selectWindow(baseTitle + ".tif" + " #" + channel_410);
	rename(baseTitle + "_410.tif");
	resetMinAndMax();
	save(dirName + getTitle);

	selectWindow(baseTitle + ".tif" + " #" + channel_470);
	rename(baseTitle + "_470.tif");
	resetMinAndMax();
	save(dirName + getTitle);

	TITLES = getList("image.titles");
}

macro "Make Masks [1]" {
	titles = getList("image.titles");
	Dialog.create("Create Masks");
	Dialog.addMessage("Select the image to create mask.");
	Dialog.addMessage("Rotations will be based on this mask.");
	Dialog.addChoice("Image", titles);
	Dialog.show();

	selectedTitle = Dialog.getChoice();
	selectWindow(selectedTitle);
	makeMask();
}

macro "Measurements [2]" {
	maskTitle = getMaskTitle();
}

macro "Test Threshold [t]" {
	threshold();
}

///////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

function getMaskTitle() {
	// First, loop through all titles and try to find one that ends with "-SEG"
	titles = getList("image.titles");
	for (i = 0; i < lengthOf(titles); i++) {
		title = titles[i];
		bName = baseName(title);
		if (endsWith(bName, "-SEG")) {
			return title;
		}
	}

	// Otherwise, we need to ask a human.
	showMessageWithCancel("No Mask Found", "No masks are open. Select one from file, or quit.")
	path = File.openDialog("Select a mask file");
	open(path);
	title = getTitle();
	return title;
}

function measureStack(title, maskTitle) {
	
}

function makeMask() {
	// Assumes that the window to run this on is open
	threshold();
	setThreshold(256, 600000);
	run("Convert to Mask...", "default");
}

function getTitles() {
	if (lengthOf(TITLES) == 0) {
		return getList("image.titles");
	} else {
		return TITLES;
	}
}

// Removes the file extension from the window name
// foobar.tif -> foobar
// foobar.txt -> foobar
function baseName(windowName) {
	baseTitle = split(windowName, '.');
	bName = baseTitle[0];
	return bName;
}

function threshold() {
	origTitle = getTitle;
	AUTOFLORESCENCE_THRESH = 1000;
	PHARYNX_THRESH = 1000;
	
	run("Duplicate...", "duplicate");
	dupID = getImageID();

	modifier = "stack";

	// First Threshold, uses edge gradients to remove autoflourescence in the gut
	run("Find Edges", modifier);
	setThreshold(AUTOFLORESCENCE_THRESH, 60000);
	setThreshold(PHARYNX_THRESH,60000);
	run("Analyze Particles...", "size=100-Infinity show=Masks exclude" + modifier);
	run("Dilate", modifier);
	run("Dilate", modifier);
	run("Dilate", modifier);
	run("Fill Holes",modifier);
	run("Erode", modifier);
	run("Erode", modifier);
	run("Erode", modifier);
	run("Divide...", "value=255 " + modifier);
	firstMaskTitle = getTitle();
	resetThreshold();
	selectWindow(origTitle);
	imageCalculator("multiply create 32-bit " + modifier, firstMaskTitle, origTitle);
	firstSegTitle = getTitle();
	
	// Second Threshold, removes background pixels left over from first mask
	setThreshold(PHARYNX_THRESH, 60000);
	run("Analyze Particles...", "size=200-Infinity pixel show=Masks " + modifier);
	run("Fill Holes", modifier);
	run("Divide...", "value=255 " + modifier);
	maskTitle = getTitle();
	imageCalculator("multiply create 32-bit " + modifier, maskTitle, origTitle);
	resetMinAndMax();
	segTitle = getTitle();
	
	// Close extraneous windows
	selectImage(dupID);
	close();
	selectWindow(firstMaskTitle);
	close();
	selectWindow(maskTitle);
	close();
	selectWindow(firstSegTitle);
	close();
	
	selectWindow(segTitle);
	rename(segName(origTitle));
	
	// Set LUT to be grayscale
	reds = newArray(256); 
	greens = newArray(256); 
	blues = newArray(256);
	for (i=0; i<256; i++) {
	    reds[i] = i;
	    greens[i] = i;
	    blues[i] = i;
	}
}

function segName(windowName) {
	bName = baseName(windowName);
	return bName + "-SEG.tif";
}
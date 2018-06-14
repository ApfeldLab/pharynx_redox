run("Set Scale...", "distance=0 global");
origTitle = getTitle();
run("Duplicate...", "duplicate");
dupID = getImageID();

// set to "stack" for stack, "slice" for individual image
modifier = getArgument();
if (lengthOf(modifier) == 0) {
	modifier = "stack";	 
}


// First Threshold, gets rid of autoflourescence
//run("Gaussian Blur...", "sigma=2");
run("Find Edges", "stack");
setThreshold(1000, 60000);
run("Analyze Particles...", "size=200-Infinity pixel show=Masks " + "stack");
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
setThreshold(1200, 60000);
run("Analyze Particles...", "size=200-Infinity pixel show=Masks " + modifier);
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
rename(origTitle + '-SEG');

// Set LUT to be grayscale
reds = newArray(256); 
greens = newArray(256); 
blues = newArray(256);
for (i=0; i<256; i++) {
    reds[i] = i;
    greens[i] = i;
    blues[i] = i;
}
//setLut(reds, greens, blues);

origTitle = getTitle();
run("Duplicate...", "duplicate");
dupID = getImageID();

// First Threshold, gets rid of autoflourescence
run("Find Edges", "stack");
setThreshold(4000, 60000);
run("Analyze Particles...", "size=200-Infinity show=Masks stack");
run("Dilate", "stack");
run("Dilate", "stack");
run("Dilate", "stack");
run("Fill Holes", "stack");
run("Erode", "stack");
run("Erode", "stack");
run("Erode", "stack");

firstMaskTitle = getTitle();
for (i=0; i <nSlices; i++) {
	setSlice(i+1);
	run("Max...", "value=1");
}
selectImage(dupID);
close();
resetThreshold();
selectWindow(origTitle);
imageCalculator("multiply create 32-bit stack", firstMaskTitle, origTitle);
firstSegTitle = getTitle();

// Second Threshold, removes background pixels left over from first mask
setThreshold(1800, 60000);
run("Analyze Particles...", "size=200-Infinity show=Masks stack");
maskTitle = getTitle();
for (i=0; i <nSlices; i++) {
	setSlice(i+1);
	run("Max...", "value=1");
}
imageCalculator("multiply create 32-bit stack", maskTitle, origTitle);
resetMinAndMax();
segTitle = getTitle();

selectWindow(firstMaskTitle);
close();
selectWindow(maskTitle);
close();
selectWindow(firstSegTitle);
close();

selectWindow(segTitle);
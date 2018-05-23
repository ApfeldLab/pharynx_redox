// Posterior -> Anterior Alignment
run("Clear Results");
run("Set Measurements...", "area shape redirect=None decimal=7");
selectWindow("PA-"+segTitle);
setThreshold(256, 60000);

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

flipIfNecessary("PA-"+origTitle);
flipIfNecessary("PA-"+segTitle);

function flipIfNecessary(title, shouldFlipArray) {
	selectWindow(title);
	for (j = 0; j < nSlices; j++) {
		setSlice(j + 1);
		if (shouldFlipArray[j]) {
			run("Flip Horizontally", "slice");
		}
	}
}
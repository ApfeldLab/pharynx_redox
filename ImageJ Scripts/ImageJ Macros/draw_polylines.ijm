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
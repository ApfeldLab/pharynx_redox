width = getWidth;
height = getHeight;

for (j = 0; j < nSlices; j++) {
	setSlice(j+1);
	
	run("Clear Results");
	snapshot;
	changeValues(0, 1200, 0);
	// Get Center of Masses
	run("Set Measurements...", "center");
	for (i = 0; i < width; i++) {
		makeRectangle(i, 0, 1, height);
		run("Measure");
	}

	firstNonNan = 0;
	for (i = 0; i < nResults; i++) {
		if (!isNaN(getResult("YM", i))) {
			firstNonNan = i;
			break;
		}
	}

	lastNonNan = 0;
		for (i = nResults - 1; i >= 0; i--) {
		if (!isNaN(getResult("YM", i))) {
			lastNonNan = i;
			break;
		}
	}

	for (i = 0; i < nResults; i++) {
		if (isNaN(getResult("YM", i))) {
			setResult("YM", i, 0);
		}
	}
	
	Xs = newArray(nResults);
	Ys = newArray(nResults);
	
	for (i = 0; i < nResults; i++) {
		Xs[i] = i;
		Ys[i] = getResult("YM", i);
	}
	
	XTrimmed = Array.slice(Xs, firstNonNan + 3, lastNonNan - 3);
	YTrimmed = Array.slice(Ys, firstNonNan + 3, lastNonNan - 3);

	Fit.doFit("3rd Degree Polynomial", XTrimmed, YTrimmed);

	NPOINTS = 5;
	XFit = newArray(NPOINTS);
	YFit = newArray(NPOINTS);
	
	for (i = 0; i < NPOINTS; i++) {
		XFit[i] = (width / NPOINTS) * i + (width / (NPOINTS) / 2);
		YFit[i] = Fit.f(XFit[i]);
	}
	
	makeSelection(6, XFit, YFit);
	run("Fit Spline");
	setLineWidth(3); 
	roiManager("add");
	reset;
}

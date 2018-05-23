title = getTitle;
for (i = 0; i < roiManager("count"); i++) {
	roiManager("select", i);
	getSelectionCoordinates(x, y);
	for (j = 0; j < lengthOf(x); j++) {
		setResult("X"+i, j, x[j]);
		setResult("Y"+i, j, y[j]); 	
	}
}
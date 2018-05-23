for (i=800; i < 900; i++) {
	setSlice(i+1);
	run("Max Inscribed Circles", "get minimum=1.5 minimum_0=0.50 closeness=5");
	
}
roiManager("sort");
idxs = newArray(roiManager("count") - 100);
for (i=0; i < roiManager("count") - nSlices; i++) {
	idxs[i] = i;
}
roiManager("select", idxs);
roiManager("delete");
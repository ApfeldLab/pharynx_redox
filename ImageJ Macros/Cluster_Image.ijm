macro "Cluster Images" {
	title = getTitle();
	path = getDirectory("image") + title;
	N_CLUSTERS = 3;
	result = exec("/Users/sean/.virtualenvs/wormAnalysis/bin/python /Users/sean/code/wormAnalysis/cluster.py " + path + " " + N_CLUSTERS);
	print(result);
	lines = split(result,"\n");
	setFont("SansSerif", 9);
	for (i = 0; i < lines.length; i++) {
		run("Make Substack...", "slices="+lines[i]);
		run("Grays");
		for (j = 0; j < nSlices; j++) {
			setSlice(j+1);
			orig_idx = split(lines[i], ",");
			drawString(orig_idx[j], 0, 10, "white");
		}
		selectWindow(title);
	}
}
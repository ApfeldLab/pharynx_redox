dirName = getInfo("image.directory");
run("Clear Results");
for (i = 0; i < roiManager("count"); i++) {
	roiManager("select", i);
	profile = getProfile();
	for (j=0; j<profile.length; j++) {
		setResult(""+i, j, profile[j]);
	}
}
updateResults;
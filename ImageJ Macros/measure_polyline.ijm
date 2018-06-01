dirName = getInfo("image.directory");
title = getTitle;
baseTitle = split(title, '.');
baseTitle = baseTitle[0];

run("Clear Results");
for (i = 0; i < roiManager("count"); i++) {
	roiManager("select", i);
	profile = getProfile();
	for (j=0; j<profile.length; j++) {
		setResult(""+i, j, profile[j]);
	}
}

updateResults;
selectWindow("Results");
save(dirName + baseTitle + "_intensities.txt");
run("Close");
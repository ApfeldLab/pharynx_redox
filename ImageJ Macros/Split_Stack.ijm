dirName = getInfo("image.directory");
baseTitle = split(getTitle, ".");
baseTitle = baseTitle[0];

Dialog.create("Split Channels");
Dialog.addMessage("Enter the channels, and comma-separate if >1 (e.g. 2,3,5)")
Dialog.addString("TL Channel(s)", "1");
Dialog.addString("470 Channel(s)", "2");
Dialog.addString("410 Channel(s)", "3");
Dialog.show();

strListTL = Dialog.getString();
strList470 = Dialog.getString();
strList410 = Dialog.getString();

listTL = split(strListTL, ",");
list470 = split(strList470, ",");
list410 = split(strList410, ",");

n_channels = lengthOf(listTL) + lengthOf(list470) + lengthOf(list410);

run("Deinterleave", "how=n_channels");

for (i=0; i < lengthOf(listTL); i++) {
	selectWindow(baseTitle + ".tif" + " #" + listTL[i]);
	rename(baseTitle + "_TL.tif");
	save(dirName + getTitle);
}

for (i=0; i < lengthOf(list470); i++) {
	selectWindow(baseTitle + ".tif" + " #" + list470[i]);
	rename(baseTitle + "_470.tif");
	save(dirName + getTitle);
}

for (i=0; i < lengthOf(list410); i++) {
	selectWindow(baseTitle + ".tif" + " #" + list410[i]);
	rename(baseTitle + "_410.tif");
	save(dirName + getTitle);
}
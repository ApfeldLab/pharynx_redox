macro "Label 0 [0]" {
	setMetadata("Label", "0");
	run("Next Slice [>]");
}
macro "Label 1 [1]" {
	setMetadata("Label", "1");
	run("Next Slice [>]");
}
macro "Label 2 [2]" {
	setMetadata("Label", "2");
	run("Next Slice [>]");
}
macro "Label 3 [3]" {
	setMetadata("Label", "3");
	run("Next Slice [>]");
}
macro "Label 4 [4]" {
	setMetadata("Label", "4");
	run("Next Slice [>]");
}

macro "Collect Labels [c]" {
	run("Clear Results");
	for (i = 0; i < nSlices; i++) {
		setSlice(i+1);
		setResult("Label", i, getMetadata("Label"));
	}
}
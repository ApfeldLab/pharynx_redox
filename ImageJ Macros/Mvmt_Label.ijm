macro "label 0 [0]" {
	run("Set Label...", "label=0");
	nextSlice();
}

macro "label 1 [1]" {
	run("Set Label...", "label=1");
	nextSlice();
}

macro "label 2 [2]" {
	run("Set Label...", "label=2");
	nextSlice();
}

macro "label 3 [3]" {
	run("Set Label...", "label=3");
	nextSlice();
}

macro "collect labels [c]" {
	for (i = 1; i <= nSlices; i++) {
		setSlice(i);
		setResult("Movement Label", i-1, getInfo("slice.label"));
	}
	updateResults();
}

function nextSlice() {
	setSlice(minOf(getSliceNumber() + 1, nSlices));
}
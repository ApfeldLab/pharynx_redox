//var movementAnnotations = newArray(nSlices);

macro "Label 0 [0]" {
	labelMovement(0);
}

macro "Label 1 [1]" {
	labelMovement(1);
}

macro "Label 2 [2]" {
	labelMovement(2);
}

macro "Label 3 [3]" {
	labelMovement(3);
}

function labelMovement(label) {
	currentSlice = getSliceNumber();
	setMetadata("Label", label);
	if (currentSlice < nSlices) {
		setSlice(currentSlice + 1);	
	}
}

macro "Save Movement Annotations To Disk [S]" {
	imageDir = getDirectory("image");
	annotationsFullPath = imageDir + File.separator + "annotations.csv";

	if (File.exists(annotationsFullPath)) {
		File.delete(annotationsFullPath);
	}

	currentSlice = getSliceNumber();
	for (i = 0; i < nSlices; i++) {
		setSlice(i + 1);
		label = getMetadata("Label");
		File.append(toString(label), annotationsFullPath);
	}
	setSlice(currentSlice);
}

macro "Print Movement Annotations To Check [P]" {
	currentSlice = getSliceNumber();
	movementAnnotations = newArray(nSlices);
	for (i = 0; i < nSlices; i++) {
		setSlice(i + 1);
		label = getMetadata("Label");
		movementAnnotations[i] = label;
	}
	setSlice(currentSlice);
	Array.show(movementAnnotations);
}

macro "Spread Windows" {
	titlesByHeight = getImageWindowTitles();
	titlesByWidth = Array.copy(titlesByHeight);
	
	sortWindowTitlesByHeight(titlesByHeight);
	sortWindowTitlesByWidth(titlesByWidth);

	
}

function getImageWindowTitles() {
	titles = newArray(nImages()); 
  	for (i=1; i<=nImages(); i++) { 
    	selectImage(i); 
    	titles[i-1] = getTitle();
    }
    return titles;
}

// Descending Order
function sortWindowTitlesByHeight(windowTitles) {
	titles = Array.copy(windowTitles);
	// Initialize list of key/value where key=title, value=height
	for (i = 0; i < titles.length; i++) {
		selectImage(titles[i]);
		height = 0;
		getLocationAndSize(None, None, None, height);
		List.set(titles[i], height);
	}

	// Insertion Sort
	i = 1;
	while (i < titles.length) {
		j = i;
		while (j > 0) {
			if (List.get(titles[j-1]) < List.get(titles[j])) {
				// Swap j-1 and j
				temp = titles[j];
				titles[j] = titles[j-1];
				titles[j-1] = temp;
				j--;
			} else {
				break;
			}
		}
		i++;
	}
}

function sortWindowTitlesByWidth(titles) {
	// Initialize list of key/value where key=title, value=height
	for (i = 0; i < titles.length; i++) {
		selectImage(titles[i]);
		width = 0;
		getLocationAndSize(None, None, width, None);
		List.set(titles[i], width);
	}


	// Insertion Sort
	i = 1;
	while (i < titles.length) {
		j = i;
		while (j > 0) {
			if (List.get(titles[j-1]) < List.get(titles[j])) {
				// Swap j-1 and j
				temp = titles[j];
				titles[j] = titles[j-1];
				titles[j-1] = temp;
				j--;
			} else {
				break;
			}
		}
		i++;
	}
}
//makeRectangle(0, 150, 110, 110); // analysis images gut
//makeRectangle(0, 21, 273, 20); // analysis images gut
// 

setThreshold (2000, 60000);
run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=�m global");
run("Set Measurements...", "area centroid center perimeter bounding fit shape feret's median skewness kurtosis scientific redirect=None decimal=7");
run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display add stack");

run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7"); 
run("Measure Stack");

macro "Measure stack [1]"{    
 run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7"); 
 setSlice(1); 
   for (i=0; i<nSlices; i++) {
	run("Measure");
	run("Next Slice [>]");
   }
beep();
}

macro "SubtractMode_PerSlice [2]"{    
 run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7"); 
 setSlice(1); 
   for (i=0; i<nSlices; i++) {
	run("Measure");
	Background = getResult("Mode");
	run("Subtract...", "value=Background");
	run("Measure");
       run("Next Slice [>]");
   }
beep();
}


macro "Sort_by_Excitation [0]" {
		frames = nSlices/5
		run("Stack to Hyperstack...", "order=xyczt(default) channels=5 slices=1 frames=frames display=Color");
		run("Split Channels");
	}


// Requires ROI manager and Results windows open

macro "Thresholded_Particles [3]" { 
 	setSlice(1); 
	run("Threshold...");
	setThreshold (1200, 60000);
	//run("Set Measurements...", "area mean standard modal min redirect=None decimal=7"); // for other tissues
	run("Set Measurements...", "area mean standard modal min centroid center bounding fit redirect=None decimal=7");
	run("Analyze Particles...", "size=200-Infinity circularity=0.00-1.00 show=Nothing display add stack");

if(roiManager("count") != nSlices){
	showMessage("Segmentation Error: ROIs count is different from Number of slices");
	}
beep();
}

macro "Roi to stack [4]"{
//if ( !(isOpen("ROI Manager"))| !(isOpen("Results")) ) {
//exit("ROI manager and Results window required")

	name = getInfo("image.filename");
	run("Set Measurements...", "area mean standard modal min centroid center bounding fit redirect=None decimal=7");       
	roiManager("Measure");
	name = getInfo("image.filename");
	type=bitDepth();
	newImage("ROI_stack", type, 100, 100, nResults);
	selectWindow(name);
	for (i=0; i<nResults; i++) {  			//paste ROIs to an empty stack
		selectWindow(name);
		Angle =  getResult('Angle', i);	
		roiManager("Select", i);
	  	run("Copy");
		selectWindow("ROI_stack");
		setSlice(i+1);
	  	run("Paste");
	  	run("Select None");
	  	run("Rotate... ", "angle=Angle grid=1 interpolation=Bilinear slice");
	 	
  	 }
selectWindow("ROI_stack");
resetMinAndMax();
roiManager("Deselect");
selectWindow("Results");
beep();
}
/}



macro "PA Alignement [5]"{
run("Set Measurements...", "area shape redirect=None decimal=7");
setSlice(13); 
setThreshold(200, 60000);
   for (i=1; i<nSlices+1; i++) {		
	//makeRectangle(0, 29, 37, 45);
	makeRectangle(0, 30, 40, 45); // for old
	run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing display");	run("Select None");	

	//makeRectangle(62, 29, 37, 45);
	makeRectangle(59, 29, 37, 45); // for old
        run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing display");	run("Select None");	

	leftEnd= getResult("Round", 0);
	rightEnd= getResult("Round", 1);	
	if (leftEnd< rightEnd){
			   run("Flip Horizontally", "slice");
		}
	selectWindow("Results");
	run("Close");
       run("Next Slice [>]");
   }
	run("Set Measurements...", "area mean standard modal min centroid center bounding fit redirect=None decimal=7");
beep();
}

//run("Flip Horizontally", "slice");



macro "7 Points [7]"{    
run("Set Measurements...", "area centroid center redirect=None decimal=7");
setSlice(1); 
setThreshold(200, 60000);
  for (i=0; i<100; i++) {	
	run("Select None");	
	makeRectangle(20, 35, 16, 32);
       run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing");

	run("Select None");	
	makeRectangle(39, 35, 3, 32);
       run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing");

	run("Select None");	
	makeRectangle(46, 35, 3, 32);
       run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing");

	run("Select None");	
	makeRectangle(51, 35, 10, 32);
       run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing");

	run("Select None");	
	//makeRectangle(62, 35, 11, 32);
	makeRectangle(67, 35, 11, 32); // for old
       run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing display");

	x1= getResult("X", 0);
	y1= getResult("Y", 0);
	x2= getResult("X", 1);
	y2= getResult("Y", 1);
	x3= getResult("X", 2);
	y3= getResult("Y", 2);
	x4= getResult("X", 3);
	y4= getResult("Y", 3);
	x5= getResult("X", 4);
	y5= getResult("Y", 4);
	
	m1=(y1-y2)/(x1-x2);
	x0=(17);
	y0=(-m1*x0 + y1);

	m2=(y5-y4)/(x5-x4);
	x6=(82);
	y6= (m2*x6 + y5);

	setLineWidth(2);
	makeLine(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6);
	setLineWidth(2);
	run("Fit Spline");
	roiManager("Add");
	selectWindow("Results");
	run("Close");
        run("Next Slice [>]");
  	 }
 run("Set Measurements...", "area mean standard modal min centroid median display redirect=None decimal=7"); 
beep();
}

macro "Measure Morphology [8]"{
setThreshold(200, 60000);
run("Set Scale...", "distance=1 known=2.58 pixel=1 unit=�m");// 1px=2.58�m for 4x4 ; 1px=1.29�m for 2x2
run("Set Measurements...", "area centroid center perimeter bounding fit shape feret's median skewness kurtosis scientific redirect=None decimal=7");
run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display add stack");
}

macro "Update ROI [9]"{
run("Create Selection");
roiManager("Update");
}


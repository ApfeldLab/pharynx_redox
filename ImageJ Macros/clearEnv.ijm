macro "Reset Environment" { 
	while (nImages>0) { 
    	selectImage(nImages); 
        close(); 
    }
    run("Clear Results");
    roiManager("deselect");
    roiManager("delete");
    run("Main Window [return]");
}
run("8-bit");
setAutoThreshold("Default dark");
//run("Threshold...");
setThreshold(20, 255);
//setThreshold(20, 255);
setOption("BlackBackground", true);
run("Convert to Mask", "method=Default background=Dark black");
run("Close");
run("Analyze Particles...", "size=200-Infinity show=Masks stack");

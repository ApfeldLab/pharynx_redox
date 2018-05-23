runMacro("/Users/sean/code/wormAnalysis/ImageJ Macros/threshold.ijm", "stack");
setThreshold (25, 60000);
run("Analyze Particles...", "size=200-Infinity pixel circularity=0.00-1 show=Masks stack");
origTitle = getTitle;
setThreshold (1200, 60000);
run("Analyze Particles...", "size=200-Infinity circularity=0.00-1.00 show=Masks stack");
run("Divide...", "value=255 stack");
maskTitle = getTitle;
imageCalculator("multiply create 32-bit stack", maskTitle, origTitle);
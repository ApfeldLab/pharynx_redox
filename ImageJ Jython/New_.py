from ij import IJ, WindowManager

IJ.setThreshold(255, 255, "Red")
IJ.run("Analyze Particles...", "size=200 exclude add stack")
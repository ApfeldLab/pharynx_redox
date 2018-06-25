# Make Measurements

from ij import IJ, WindowManager
from ij.gui import GenericDialog, Roi
from ij.measure import ResultsTable, Measurements, Calibration
from ij.process import ImageStatistics, ImageProcessor as IP
from ij.plugin.frame import RoiManager
from ij.plugin.filter import ParticleAnalyzer as PA, PlugInFilter as PF, Analyzer

#### Helper Functions
def getMedians(imPlus):
	medians = []
	stk = imPlus.getStack()
	for i in range(stk.getSize()):
		img_stats = ImageStatistics.getStatistics(stk.getProcessor(i+1), Measurements.MEDIAN, Calibration())
		medians.append(img_stats.median)
	return medians

imTitles = WindowManager.getImageTitles()
nTitles = len(imTitles)

## DIALOG
d = GenericDialog("Measure")
d.addChoice("Mask", imTitles, "-")
d.addChoice("470 image", imTitles, "-")
d.addChoice("410 image", imTitles, "-")
d.showDialog()

imgMask = WindowManager.getImage(d.getNextChoice())
img410 = WindowManager.getImage(d.getNextChoice())
img470 = WindowManager.getImage(d.getNextChoice())

##################################################
# MEASUREMENTS
table = ResultsTable()
roiManager = RoiManager(True) # Boolean => Don't Display
PA.setRoiManager(roiManager)

# Measure Median 410 and 470
med410 = getMedians(img410)
med470 = getMedians(img470)

# Make the ROIs based on the mask
IJ.setThreshold(imgMask, 255, 255, "No Update")
IJ.run(imgMask, "Analyze Particles...", "size=200 exclude add stack") # TODO: why is there a border on the mask? This is why I am using exclude
roiManager.runCommand("Show None")

maskStk = imgMask.getStack()
im410Stk = img410.getStack()
im470Stk = img470.getStack()

im410_means = []
im470_means = []
areas = []
angles = []
xs = []
ys = []

ra = roiManager.getRoisAsArray()
for i in range(maskStk.getSize()):
	ip = im410Stk.getProcessor(i+1)
	ip.setRoi(ra[i])
	istats = ip.getStatistics()
	im410_means.append(istats.mean)

	ip = im470Stk.getProcessor(i+1)
	ip.setRoi(ra[i])
	istats = ip.getStatistics()
	im470_means.append(istats.mean)

	# TODO: Figure out units for Area
	areas.append(istats.area)
	angles.append(istats.angle)
	xs.append(istats.xCenterOfMass)
	ys.append(istats.yCenterOfMass)


###############################################
# Orient & Align
img410PA = img410.duplicate()
img470PA = img470.duplicate()
maskPA = imgMask.duplicate()

def translate(ip, x_center_obj, y_center_obj):
	x_center_im = ip.getWidth() / 2
	y_center_im = ip.getHeight() / 2
	x_offset = x_center_im - x_center_obj
	y_offset = y_center_im - y_center_obj
	ip.translate(x_offset, y_offset)

def rotateStack(imPlus, angles, interp_method):
	""" Rotate a stack of images

	`angles` is an array of angles in radians # TODO: radians or degrees?
	`interp_method` is from the ImageProcessor class (ImageProcessor.BILINEAR, etc.)
	"""
	stk = imPlus.getStack()
	for i in range(stk.getSize()):
		ip = stk.getProcessor(i+1)
		
		ip.setInterpolationMethod(interp_method)
		ip.setBackgroundValue(0)
		translate(ip, xs[i], ys[i])
		ip.rotate(angles[i])

def getFlipArray(maskPAImgPlus):
	maskPAImgPlus.show()
	rt = ResultsTable()
	PA.setResultsTable(rt)
	IJ.selectWindow(maskPAImgPlus.getTitle())
	IJ.run("Set Measurements...", "shape")
	IJ.setThreshold(maskPAImgPlus, 1, 1000, "Red")
	stk = maskPAImgPlus.getStack()
	IJ.makeRectangle(maskPAImgPlus.getWidth()/2-5,0,maskPAImgPlus.getWidth(),maskPAImgPlus.getHeight())
	
	for i in range(stk.getSize()):
		IJ.setSlice(i+1)
		IJ.run(imgMask, "Analyze Particles...", "exclude")
	rt.show("Aligning")


rotateStack(img410PA, angles, IP.BILINEAR)
rotateStack(maskPA, angles, IP.NONE)
img410PA.show()
getFlipArray(maskPA)
maskPA.show()


# Measure morphology
#TODO
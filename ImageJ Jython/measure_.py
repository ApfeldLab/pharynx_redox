from __future__ import division

import csv
import os
from collections import defaultdict
import json

from ij import IJ, WindowManager
from ij.gui import GenericDialog, Roi
from ij.measure import ResultsTable, Measurements, Calibration
from ij.process import ImageStatistics, ImageProcessor as IP, StackProcessor
from ij.plugin.frame import RoiManager
from ij.plugin.filter import ParticleAnalyzer as PA, PlugInFilter as PF, Analyzer


##################################################
# UTILITY HELPER FUNCTIONS

def getMedians(imPlus):
	medians = []
	stk = imPlus.getStack()
	for i in range(stk.getSize()):
		img_stats = ImageStatistics.getStatistics(stk.getProcessor(i+1), Measurements.MEDIAN, Calibration())
		medians.append(img_stats.median)
	return medians

def addValues(table, column, values):
	"""Add values (list) to column of table"""
	for i,value in enumerate(values):
		table.setValue(column, i, value)

def pa_title(imageTitle):
	"""Give the PA title for the given image title (including the .tif)"""
	baseTitle = imageTitle.split(".")[0]
	return baseTitle + "-PA" + ".tif"

def mask_title(imageTitle):
	"""Give the mask title for the given image title (including the .tif)"""
	baseTitle = imageTitle.split(".")[0]
	return baseTitle + "-MASK" + ".tif"

def get_processors(image_plus):
	stk = image_plus.getStack()
	for i in range(stk.getSize()):
		yield stk.getProcessor(i+1)

def applyToStack(im_plus, fun, kwargs):
	"""Apply a function to each ImageProcessor in the given ImagePlus

	The function should take in an ImageProcessor as its first argument. Other arguments can
	be passed with kwargs.
	
	e.g.
	def addN(imProc, n):
		imProc.add(n)
	
	applyToStack(an_imPlus, addN, 10)
	"""
	map(lambda im_proc: fun(im_proc, **kwargs), get_processors(im_plus))

###################################################
# SETUP

imTitles = WindowManager.getImageTitles()
nTitles = len(imTitles)

## DIALOG
d = GenericDialog("Measure")
d.addChoice("Mask", imTitles, "-")
d.addChoice("470 image", imTitles, "-")
d.addChoice("410 image", imTitles, "-")
d.addRadioButtonGroup("Binning", ["4x4", "2x2"], 1, 2, "4x4")
d.showDialog()

imgMask = WindowManager.getImage(d.getNextChoice())
img410 = WindowManager.getImage(d.getNextChoice())
img470 = WindowManager.getImage(d.getNextChoice())
binning = d.getNextRadioButton()

PARENT_DIR = img410.getOriginalFileInfo().directory

# setup data collection object
dataTable = ResultsTable()

##################################################
# INITIAL MEASUREMENTS
table = ResultsTable()
roiManager = RoiManager(True) # Boolean => Don't Display
PA.setRoiManager(roiManager)

# Measure Median 410 and 470
addValues(dataTable, 'Median 410', getMedians(img410))
addValues(dataTable, 'Median 470', getMedians(img470))

# Make the ROIs based on the mask
IJ.setThreshold(imgMask, 255, 255, "No Update")
IJ.run(imgMask, "Analyze Particles...", "size=200 exclude add stack")
roiManager.runCommand("Show None")

maskStk = imgMask.getStack()
im410Stk = img410.getStack()
im470Stk = img470.getStack()

areas = []
angles = []
xs = []
ys = []

# TODO: subtract median here?
ra = roiManager.getRoisAsArray()
for i in range(maskStk.getSize()):
	ip = im410Stk.getProcessor(i+1)
	ip.setRoi(ra[i])
	istats = ip.getStatistics()
	dataTable.setValue('Intensity410_wholePharynx', i, istats.mean)

	ip = im470Stk.getProcessor(i+1)
	ip.setRoi(ra[i])
	istats = ip.getStatistics()
	dataTable.setValue('Intensity470_wholePharynx', i, istats.mean)

	# TODO: Figure out units for Area
	dataTable.setValue('Area (Px)', i, istats.area)

	# we don't need to keep track of these in our `data` object, just need them to do the rotations
	angles.append(istats.angle)
	xs.append(istats.xCenterOfMass)
	ys.append(istats.yCenterOfMass)

roiManager.runCommand("reset")


###############################################
# Orient & Align

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
	for i, ip in enumerate(get_processors(imPlus)):
		ip.setInterpolationMethod(interp_method)
		ip.setBackgroundValue(0)
		translate(ip, xs[i], ys[i])
		ip.rotate(angles[i])
		

def getFlipArray(maskPAImgPlus):
	"""Given the PA-aligned mask image, return a list of booleans, indicating whether or not we should flip the corresponding index"""
	maskPAImgPlus.show()
	IJ.run(maskPAImgPlus, "Shape Smoothing", "relative_proportion_fds=15 absolute_number_fds=2 keep=[Relative_proportion of FDs] stack");
	
	IJ.run("Set Measurements...", "shape")
	IJ.setThreshold(maskPAImgPlus, 1, 1000, "No Update")
	stk = maskPAImgPlus.getStack()

	tableLeft = ResultsTable()
	PA.setResultsTable(tableLeft)
	IJ.makeRectangle(0,0,maskPAImgPlus.getWidth()/2-5,maskPAImgPlus.getHeight())
	IJ.run(maskPAImgPlus, "Analyze Particles...", "stack")

	tableRight = ResultsTable()
	PA.setResultsTable(tableRight)
	IJ.makeRectangle(maskPAImgPlus.getWidth()/2-5,0,maskPAImgPlus.getWidth()/2+5,maskPAImgPlus.getHeight())
	IJ.run(maskPAImgPlus, "Analyze Particles...", "stack")
	maskPAImgPlus.hide()

	IJ.resetThreshold(maskPAImgPlus)

	leftCirc = [tableLeft.getValue("Circ.", i) for i in range(stk.getSize())]
	rightCirc = [tableRight.getValue("Circ.", i) for i in range(stk.getSize())]
	return [leftCirc[i] < rightCirc[i] for i in range(len(leftCirc))]

def flipIfNecessary(flipArray, imPluses):
	for imP in imPluses:
		stk = imP.getStack()
		for i in range(stk.getSize()):
			proc = stk.getProcessor(i+1)
			if flipArray[i]:
				proc.flipHorizontal()

img410PA = img410.duplicate()
img470PA = img470.duplicate()
maskPA = imgMask.duplicate()

for imPlus in [img410PA, img470PA, maskPA]:
	imPlus.setTitle(pa_title(imPlus.getTitle().lstrip("DUP_")))

rotateStack(img410PA, angles, IP.BILINEAR)
rotateStack(img470PA, angles, IP.BILINEAR)
rotateStack(maskPA, angles, IP.NONE)
flipIfNecessary(getFlipArray(maskPA), [img410PA, img470PA, maskPA])

img410PA.show()
img470PA.show()

# MORPHOLOGICAL MEASUREMENTS
# These are done on the rotated image mask

scales = {
	"4x4": 2.58,
	"2x2": 1.29
}

morphTable = ResultsTable()
PA.setResultsTable(morphTable)

maskPA.show()
IJ.run(maskPA, "Select None", "")
IJ.run(maskPA, "Set Scale...", "distance=1 known=" + str(scales[binning]) + " pixel=1 unit=μm"); #TODO scale?
IJ.run(maskPA, "Set Measurements...", "area centroid center perimeter bounding fit shape feret's median skewness kurtosis scientific redirect=None decimal=7");
IJ.run(maskPA, "Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing stack");

morph_headings = morphTable.getColumnHeadings().strip(' ').split('\t')
print(morph_headings)
for m_heading in morph_headings[1:]: # because of formatting, the first is empty
	addValues(dataTable, m_heading, [morphTable.getValue(m_heading, i) for i in range(morphTable.size())])

def cropStack(imPlus, width=120, height=50):
	x = int(imPlus.getWidth() / 2 - width / 2)
	y = int(imPlus.getHeight() / 2 - height / 2)
	sp = StackProcessor(imPlus.getStack())
	imPlus.setStack(sp.crop(x, y, int(width), int(height)))

for imPlus in [img410PA, img470PA, maskPA]:
	cropStack(imPlus)
	IJ.saveAsTiff(imPlus, os.path.join(PARENT_DIR, imPlus.getTitle()))

dataTable.show("Measurements")
dataTable.save(os.path.join(PARENT_DIR, 'measurements.csv'))
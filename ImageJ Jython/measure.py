from __future__ import division

## Use this to recompile Jython modules to class files.
#from sys import modules
#modules.clear()

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

from pharynx_analysis.utils import *

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
img410PA = img410.duplicate()
img470PA = img470.duplicate()
maskPA = imgMask.duplicate()

for imPlus in [img410PA, img470PA, maskPA]:
	imPlus.setTitle(pa_title(imPlus.getTitle().lstrip("DUP_")))

rotateStack(img410PA, angles, xs, ys, IP.BILINEAR)
rotateStack(img470PA, angles, xs, ys, IP.BILINEAR)
rotateStack(maskPA, angles, xs, ys, IP.NONE)
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
for m_heading in morph_headings[1:]: # because of formatting, the first is empty
	addValues(dataTable, m_heading, [morphTable.getValue(m_heading, i) for i in range(morphTable.size())])

for imPlus in [img410PA, img470PA, maskPA]:
	cropStack(imPlus)
	IJ.saveAsTiff(imPlus, os.path.join(PARENT_DIR, imPlus.getTitle()))

dataTable.show("Measurements")
dataTable.save(os.path.join(PARENT_DIR, 'measurements.csv'))
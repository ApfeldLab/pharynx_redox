from ij import IJ, WindowManager
from ij.gui import GenericDialog, Roi
from ij.measure import ResultsTable, Measurements, Calibration
from ij.process import ImageStatistics, ImageProcessor as IP, StackProcessor
from ij.plugin.frame import RoiManager
from ij.plugin.filter import ParticleAnalyzer as PA, PlugInFilter as PF, Analyzer

def getMedians(imPlus):
	medians = []
	stk = imPlus.getStack()
	for i in range(stk.getSize()):
		img_stats = ImageStatistics.getStatistics(stk.getProcessor(i+1), Measurements.MEDIAN, Calibration())
		medians.append(img_stats.median)
	print('GOT MEDIANS')
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


def cropStack(imPlus, width=120, height=50):
	x = int(imPlus.getWidth() / 2 - width / 2)
	y = int(imPlus.getHeight() / 2 - height / 2)
	sp = StackProcessor(imPlus.getStack())
	imPlus.setStack(sp.crop(x, y, int(width), int(height)))

def rotateStack(imPlus, angles, xCenters, yCenters, interp_method):
	""" Rotate a stack of images

	`angles` is an array of angles in radians # TODO: radians or degrees?
	`interp_method` is from the ImageProcessor class (ImageProcessor.BILINEAR, etc.)
	"""
	stk = imPlus.getStack()
	for i, ip in enumerate(get_processors(imPlus)):
		ip.setInterpolationMethod(interp_method)
		ip.setBackgroundValue(0)
		translate(ip, xCenters[i], yCenters[i])
		ip.rotate(angles[i])
		

def getFlipArray(maskPAImgPlus):
	"""Given the PA-aligned mask image, return a list of booleans, indicating whether or not we should flip the corresponding index"""
	maskPAImgPlus.show()
#	IJ.run(maskPAImgPlus, "Shape Smoothing", "relative_proportion_fds=15 absolute_number_fds=2 keep=[Relative_proportion of FDs] stack");
	
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

def translate(ip, x_center_obj, y_center_obj):
	x_center_im = ip.getWidth() / 2
	y_center_im = ip.getHeight() / 2
	x_offset = x_center_im - x_center_obj
	y_offset = y_center_im - y_center_obj
	ip.translate(x_offset, y_offset)
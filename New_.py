from ij import IJ, ImageStack, Prefs, ImagePlus
from ij.process import ByteProcessor, ImageProcessor, StackConverter
from ij.plugin import Duplicator
from ij.plugin.filter import Binary, ParticleAnalyzer as PA

binner = Binary()
def thresh(proc):
	Prefs.blackBackground = True
	proc.findEdges()
	proc.threshold(20)

	binner = Binary()
	binner.setup("fill", None)
	binner.run(proc)
	binner.setup("dilate", None)
	binner.run(proc)
	binner.setup("erode", None)
	binner.run(proc)
	return proc

im1 = IJ.getImage()
im2 = Duplicator().run(im1);

StackConverter(im2).convertToGray8()
s = im2.getStack()


minSize = 200
maxSize = 10000
options = PA.SHOW_MASKS
p = PA(options, 0, None, minSize, maxSize)
p.setHideOutputImage(True)
stk = ImageStack(im2.getWidth(), im2.getHeight())

for i in range(s.getSize()):
	im2.setSliceWithoutUpdate(i + 1)
	ip = s.getProcessor(i+1)
	s.setProcessor(thresh(ip), i+1)

	p.analyze(im2)
	mmap = p.getOutputImage()
	binner.setup("close", None)
	binner.run(mmap.getProcessor())
	stk.addSlice(mmap.getProcessor())

impThresholded = ImagePlus("tt", stk)
IJ.run(impThresholded, "Convert to Mask", "method=default background=default black")
impThresholded.show()
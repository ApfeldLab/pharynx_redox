from ij import IJ, ImageStack, Prefs, ImagePlus
from ij.process import ByteProcessor, ImageProcessor, StackConverter
from ij.plugin import Duplicator
from ij.plugin.filter import Binary, ParticleAnalyzer as PA

import os

################
# Start Script #
################
bin = Binary()
im1 = IJ.getImage()
dup_im = Duplicator().run(im1)

#StackConverter(dup_im).convertToGray8()
dup_stack = dup_im.getStack()

# First, we threshold based on edges
IJ.setThreshold(dup_im, 2500, 100000, "Red")
for i in range(dup_im.getImageStackSize()):
	dup_stack.getProcessor(i+1).findEdges()

IJ.run(dup_im, "Make Binary", "method=Default background=Default black")
IJ.run(dup_im, "Dilate", "stack")
IJ.run(dup_im, "Fill Holes", "stack")
IJ.run(dup_im, "Erode", "stack")
IJ.run(dup_im, "Erode", "stack")

# Then, remove small particles
stk = ImageStack(dup_im.getWidth(), dup_im.getHeight())
p = PA(PA.SHOW_MASKS, 0, None, 200, 100000)
p.setHideOutputImage(True)
for i in range(dup_stack.getSize()):
	dup_im.setSliceWithoutUpdate(i + 1)
	p.analyze(dup_im)
	mmap = p.getOutputImage()
	stk.addSlice(mmap.getProcessor())

dup_im.setStack(stk)
dup_im.setSliceWithoutUpdate(1)
dup_im.show()

# Finally, Save
dup_im.setTitle(im1.getTitle().split(".")[0] + "-MASK" + ".tif")

IJ.save(dup_im, os.path.join(im1.getOriginalFileInfo().directory, dup_im.getTitle()))
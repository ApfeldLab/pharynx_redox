from ij import IJ, ImageStack, Prefs, ImagePlus
from ij.process import ByteProcessor, ImageProcessor, StackConverter
from ij.plugin import Duplicator
from ij.plugin.filter import ParticleAnalyzer as PA

import os


def mask_title(imageTitle):
    """Give the mask title for the given image title (including the .tif)"""
    baseTitle = imageTitle.split(".")[0]
    return baseTitle + "-MASK" + ".tif"


def threshold(imPlus, edgeThreshold=2500):
    mask = Duplicator().run(imPlus)
    mask_stk = mask.getStack()

    # First, we threshold based on edges
    IJ.setThreshold(mask, edgeThreshold, 100000, "No Update")
    for i in range(mask.getImageStackSize()):
        mask_stk.getProcessor(i + 1).findEdges()
    IJ.run(mask, "Make Binary", "method=Default background=Default black")

    # Now, we need to clean up the binary images morphologically
    IJ.run(mask, "Dilate", "stack")
    IJ.run(mask, "Fill Holes", "stack")
    IJ.run(mask, "Erode", "stack")
    IJ.run(mask, "Erode", "stack")

    # Finally, remove the small particles
    stk = ImageStack(mask.getWidth(), mask.getHeight())
    p = PA(PA.SHOW_MASKS, 0, None, 200, 100000)
    p.setHideOutputImage(True)
    for i in range(mask_stk.getSize()):
        mask.setSliceWithoutUpdate(i + 1)
        p.analyze(mask)
        mmap = p.getOutputImage()
        stk.addSlice(mmap.getProcessor())

    mask.setStack(stk)
    mask.setSliceWithoutUpdate(1)
    mask.setTitle(mask_title(imPlus.getTitle()))
    mask.show()
    return mask


if __name__ in ["__builtin__", "__main__"]:
    threshold(IJ.getImage(), 2000)

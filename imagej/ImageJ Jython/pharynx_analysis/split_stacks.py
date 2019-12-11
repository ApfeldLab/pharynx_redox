import sys
from ij.gui import GenericDialog
from ij import IJ, WindowManager


# If there are duplicate channels, assign numbers to all channels
# and return a list of "numbered" channels. It preserves the order of
# the original list
#
# e.g. ["TL", "410", "470", "470"] -> ["TL_1", "410_1", "470_1", "470_2"]
#      ["TL", "410", "470"] -> ["TL", "410", "470"]
# 	   ["TL", 410, "470", "410", "470"] -> ["TL_1", "410_1", "470_1", "410_2", "470_2"]
def number_duplicate_channels(channels):
    if len(channels) == len(set(channels)):
        return channels
    count = {}
    numbered = []
    for channel in channels:
        if channel in count:
            count[channel] += 1
        else:
            count[channel] = 1
        numbered.append(channel + "_" + str(count[channel]))
    return numbered


# Split the given ImagePlus
def split_and_save(imp):
    base_title = imp.getShortTitle()
    base_dir = imp.getOriginalFileInfo().directory

    d = GenericDialog("Split Stacks")
    d.addMessage("Indicate the channels in your image stack, separated by commas")
    d.addMessage("Add the channels *in the order that you took the images*")
    d.addMessage('For example, "TL, 410, 470, 410, 470"')
    d.addStringField("Channels", "TL, 470, 410", 40)
    d.showDialog()
    # exit if cancelled
    if d.wasCanceled():
        return None

    channels = number_duplicate_channels(
        [x.strip() for x in d.getNextString().split(",")]
    )

    # Check that the number of images in stack is divisible by the number of channels
    # If not, quit. Else keep going
    if imp.getNSlices() % len(channels) != 0:
        IJ.error(
            "Invalid Channel Specification",
            "The number of slices (%d) is not divisible by the number of channels (%d). Exiting."
            % (imp.getNSlices(), len(channels)),
        )
        return None
    imp.show()
    IJ.run("Deinterleave", "how=%d" % len(channels))

    for i, img_id in enumerate(WindowManager.getIDList()):
        channel_i = WindowManager.getImage(img_id)
        channel_i.setTitle(base_title + "_" + channels[i])
        IJ.saveAs(channel_i, "tif", base_dir + channel_i.getTitle())
    IJ.showMessage("Saved image stacks as separate channels in %s" % base_dir)


# implus = IJ.openImage("/Users/sean/Desktop/2018_05_24/eat5_2do_05_24.tif")
# split_and_save(implus)

split_and_save(IJ.getImage())

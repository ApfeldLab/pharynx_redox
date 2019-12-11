import os
import csv

from ij.gui import GenericDialog, YesNoCancelDialog
from ij.io import OpenDialog


def get_metadata_table_filepath():
    """Get a the metadata template filename"""
    d = GenericDialog("Metadata Template")
    d.addMessage("Please choose a metadata template")
    d.enableYesNoCancel("OK", "Cancel")
    d.hideCancelButton()
    d.showDialog()
    if d.wasOKed():
        od = OpenDialog("Metadata File")
        return os.path.join(od.getDirectory(), od.getFileName())
    else:
        return


def get_metadata_rt():
    """Get the metadata template ResultsTable"""
    fn = get_metadata_table_filepath()
    with open(fn, "r") as f:
        colstr = f.read()
    cols = colstr.split(",")


def main():
    get_metadata_rt()


if __name__ in ["__builtin__", "__main__"]:
    main()

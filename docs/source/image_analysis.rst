==============
Image Analysis
==============

This document goes over how to perform the ratiometric flourescence image analysis.

Environment & Installation
##########################

Software Required:
    * Git
    * FIJI
    * Python3
    * MATLAB
    * JMP

MacOS
#####

Before we install other stuff, we'll need to install `Homebrew <https://brew.sh/>`_. It's a package manager for MacOS. Essentially, it makes installing many programs as easy as::

    $ brew install python3

Git
***
To check if git is installed, run::

    $ which git

if you get output that looks like::

    /usr/bin/git

You're good to go. Otherwise, run::

    $ brew install git

FIJI
****

Download fiji here: https://downloads.imagej.net/fiji/latest/fiji-macosx.zip

Extract the zip (it might auto-extract), then drag the Fiji.app file to your /Applications directory.

Once it's installed, install these plugins (`how to install plugins <>`_)
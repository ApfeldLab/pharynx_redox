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

Once it's installed, install these plugins (`how to install plugins <>`_):
    
    * https://imagej.net/Shape_Smoothing


MATLAB
******

`Follow Northeastern's instructions <https://northeastern.service-now.com/kb_view.do?sysparm_article=KB0012561>`_.

JMP
***

You can download JMP through the myNortheastern portal.


Downloading The Code
####################

Open Terminal.app. Change to a directory that you would like to download your code in. I suggest making a
directory in your home folder called code, and downloading the analysis there. The process looks like this::

    $ mkdir -p ~/code
    $ cd ~/code

Once you're in whatever folder you'd like to download the code in, run the following command (you need to have git installed for this)::

    $ git clone https://github.com/half-adder/wormAnalysis.git
    $ cd wormAnalysis

Now your terminal should be pointed to the wormAnalysis directory with all the code in it. To verify, run::

    $ open .

This should open the directory and show you all of the nice files and folders.

Running the Analysis
####################

Setup
*****

You should have a single .tif file from Metamorph. Place that .tif in a new directory in the shared folder, SharedFolder > Image Analysis > data > EXPERIMENT_NAME > EXPERIMENT_NAME.tif.

.. warning::
    Make sure that the .tif file is named *exactly* the same as the folder it resides in.

ImageJ Processing
*****************

Open your .tif in ImageJ 

.. note::
    I've noticed that opening directly from the shared directory makes things slow because it has to transfer data over the network. Instead, what I like to do is copy the EXPERIMENT_NAME directory to my desktop, then copy it back to the shared directory when I'm done with the analysis).


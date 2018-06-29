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
=====

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

Open ``Terminal.app``. Change to a directory that you would like to download your code in. I suggest making a
directory in your home folder called code, and downloading the analysis there. The process looks like this::

    $ mkdir -p ~/code
    $ cd ~/code

Once you're in whatever folder you'd like to download the code in, run the following command (you need to have git installed for this)::

    $ git clone https://github.com/half-adder/wormAnalysis.git
    $ cd wormAnalysis

Now your terminal should be pointed to the ``wormAnalysis`` directory with all the code in it. To verify, run::

    $ open .

This should open the directory and show you all of the nice files and folders.

ImageJ Processing
#################

Setup
*****

You should have a single .tif file from Metamorph. Place that .tif in a new directory in the shared folder, ``SharedFolder > Image Analysis > data > EXPERIMENT_NAME > EXPERIMENT_NAME.tif``.

.. warning::
    Make sure that the .tif file is named *exactly* the same as the folder it resides in. The MATLAB scripts expect this to be the case, and they throw up if it isn't. Then you'll have to go through and manually rename a bunch of files, so it's easier to just make sure it's right at this early stage.

Open your .tif in ImageJ 

.. note::
    I've noticed that opening directly from the shared directory makes things slow because it has to transfer data over the network. Instead, what I like to do is copy the ``EXPERIMENT_NAME`` directory to my desktop, then copy it back to the shared directory when I'm done with the analysis).

Now we'll install the plugins. #TODO

Split Stack
===========

The first step in our processing is to split the stack into different channels. Run the ``Split Stack`` command (mapped to ``0``).

Thresholding
============
Once the stacks are split, we need to make a threshold. The thresholding is important as the mask we generate in this step allows us to make the rotations and draw the polyline. So it's important we get a good mask here! Run the ``Threshold`` command (mapped to ``1``).

In my experience the thresholding is pretty good, but we need it to be really good, so it's important that you check that:
    
    1. There are the right number of slices
    2. Each slice contains an object that looks like a pharynx

Sometimes, the isthmus of the pharynx is left out of the mask. If this happens, you can just paint it in with ImageJ's paintbrush tool. The foreground is set to ``255``.

Othertimes, the threshold picks up some autoflourescence in the gut. If this happens, just paint over the errant foreground with the paintbrush set to ``0``.

Finally, if the thresholding macro just utterly fails, just resort to a static threshold. Run ``Image > Adjust > Threshold`` and set the slider so that you get pharynx-looking objects. Then run ``Process > Binary > Convert To Mask``, and do your painting on the resulting images.

Once you have a solid mask, you'll need to save it. If you ran the ``Threshold`` macro, the image title is correct and you just need to save it as is. However, if you used the alternative method, you'll need to specify the name manually. The name doesn't *really* matter, but I like to use the following convention: ``EXPERIMENT_NAME_CHANNEL-MASK.tif``. Save it in your experiment directory.


Pharynx-Level Measurements
==========================

This part's easy! Just run the ``Measure`` macro. This does a few things:

    1. Measure and subtract the Median from each channel.
    2. Measure the mean intensity inside the pharynx of each channel.
    3. Orient and align the pharynx.
    4. Measure pharyngeal morphology.

It will automagically save the measurements and oriented images to the experiment directory. The measurements are saved to ``EXPERIMENT_NAME_measurements.csv``. The oriented and aligned images are saved to ``EXPERIMENT_NAME_CHANNEL-PA.tif``.

Polyline Measurements
=====================

We find the medial axis of the pharynx with the oriented and aligned mask of the pharynx.

    1. Run the ``Draw Polyline`` macro
    2. Select the PA channel you would like to measure
    3. Click through the ROIs, they should be updating on the selected PA channel
    4. Adjust the polylines as necessary
    5. Run the ``Save Polyline`` macro (this saves both the coordinates and intensities to the ``EXPERIMENT_NAME`` directory)

Coalescing ImageJ data in JMP
#############################

Open the ``wormAnalysis > scripts > coalesce.jsl`` script in JMP. This script sucks in all the data that we just generated with ImageJ and, as the name implies, coalesces it into a single data table.

    1. Click ``Run Script`` in the top-left
    2. Select your ``EXPERIMENT_NAME`` directory that we just did all of our analysis in.
=======================
Developer Documentation
=======================

Dev Environment
###############


Text Editors
************

Of course, if you'd like to write code, you need something to write code *in*. If you
don't already have a preffered code editor, I would recommend downloading
`Pycharm <https://www.jetbrains.com/pycharm/>`_ (the free version works just fine, but
you can also get a free "professional" license by following the directions
`here <https://www.jetbrains.com/community/education/#students>`_. PyCharm is nice
because it has built-in auto-complete while you're typing, and has a very nice
debugger, which is critical for code development.

Configuring Code Formatting
---------------------------

This project uses the `Black <https://black.readthedocs.io/en/stable/index.html>`_
package to format all code automatically. This is so that all of our code "looks" the
same. We'll set up your editor so that it formats the file using Black every time you
save. Follow the instructions `here <https://black.readthedocs
.io/en/stable/editor_integration.html#editor-integration>`_ to configure ``black`` to
work with your text editor.

Configure Pycharm Project Interpreter
-------------------------------------

In order to do it's fancy auto-complete and other features, PyCharm needs to know
which Python environment you will be using. Follow their directions `here
<https://www.jetbrains.com/help/pycharm/configuring-python-interpreter.html>`_ to set
this up. If you are using virtual environments, use the location of that virtual
environment for this step.

This section is meant to introduce new developers to the architecture / design of this
system.

Data Structures
===============
There are a few key data structures with which one needs to be familiar in order
to work on this code. We deal with four main types of data. Those are:

    - images
    - intensity profile data in two forms:

        - raw, extracted from the images
        - functional representations of the raw data
    - tabular "summary" data

Each of these types of data are naturally internally represented in different ways.
We will go over each individually.

Images
------
Images are represented internally as numerical matrices, specifically as numpy
matrices. Thinking of a single black-and-white image, each pixel is represented by a
single number. The range of that number depends on the data type used to store the
image internally. In our case, the raw images from our microscope come to us as
16-bit unsigned integers. This means that each pixel can take any value in the range
[0, 2^16].

Still thinking of a single black-and-white image, we can access individual pixels by
*indexing* into the rows and columns of that image. For example, to get the value 10
pixels down and 3 pixels across (starting from the top left), we would use the
following notation::

    img[10, 3]

This graphic shows examples of more advanced matrix indexing.

.. image:: _static/numpy_indexing.png

Many images can be 'stacked' on top of each other, as if they were sheets of paper.
Numpy handles this case as well. All that needs to be done is to add another
dimension to our matrices. Now they are three-dimensional, with the first dimension
indicating which sheet of paper we are dealing with, and the second and third
indicating the rows and columns as before.



.. image:: _static/numpy_3d.png

For an even deeper dive on indexing, see the `numpy indexing documentation
<https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html>`_.

Writing New Code
================

Formatting
----------
The python in this code-base is formatted via the `Black <https://black.readthedocs
.io/en/stable/>`_ package. From their docs::

    By using Black, you agree to cede control over minutiae of hand-formatting. In
    return, Black gives you speed, determinism, and freedom from pycodestyle nagging
    about formatting. You will save time and mental energy for more important matters.

    Black makes code review faster by producing the smallest diffs possible. Blackened
    code looks the same regardless of the project you’re reading. Formatting becomes
    transparent after a while and you can focus on the content instead.

Please review their `documentation <https://black.readthedocs
.io/en/stable/editor_integration.html>`_ to set up your IDE to auto-format your code
with Black.

Documentation
-------------
All docstrings should be formatted in the `Numpy docstrings format <https://numpydoc
.readthedocs.io/en/latest/format.html>`_.

Adding Packages
===============
Package management is orchestrated through Anaconda_. To install a new package, use::

    $ conda add <package>

To update the list of required packages, use::

    $ conda list --explicit > conda-spec-file.txt

Building Documentation
======================
This documentation is written in RST files, and built using Sphinx. All documentation
should be written in ``docs/source``. This documentation is auto-built and uploaded to 
readthedocs on push. To build the documentation as HTML on your local machine, use::

    $ cd docs
    $ make html

The output is then in ``docs/build/html``

.. _Anaconda: https://docs.conda.io/projects/conda/en/latest/index.html
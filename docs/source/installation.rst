.. _installation:

############
Installation
############

Basic Installation
==================

Ensure that you have python installed on your system.::

    $ python --version
    Python 3.7.6

.. warning::
    Pharedox is not tested for Python < 3.7


Before moving forward, I would suggest that you use a virtual environments. Virtual
environments keep packages that PhaRedox needs separate from your system's global
python environment. See the instructions `here <https://python-guide-cn.readthedocs
.io/en/latest/dev/virtualenvs.html>`_. Note that this is *not required*, but is
helpful in the long run.

.. note::
    Make sure to call your virtual environment something meaningful, like ``pharedox``.

If you did set up a virtual environment, make sure to activate it whenever you want
to use PhaRedox. The process should look something like this::

    $ workon pharedox
    (pharedox) $

If you think that you will want to edit the source code, follow the instructions on
`How to set up a development environment`_. Otherwise, after you have set up a
virtual environment (or if you skipped that), simply execute the following line::

    pip install pharedox

MATLAB
======

If you would like to use the `1D Profile Registration` features, you will need to
set up the MATLAB python engine. Ensure that MATLAB is installed before continuing
(your university probably has instructions on how to do this).

.. warning::
    PhaRedox was developed using ``MATLAB_R2019a``. Other versions are not guaranteed
    to work, though you are free to try.


Add PhaRedox files to MATLAB path
*********************************
Open up MATLAB. Look for ``Set Path`` in the ``Home`` tab. Click it, then in the dialog,
click ``Add with subfolders``, and navigate to the PhaRedox source directory and select
the ``matlab`` folder.

Install MATLAB engine in Python
*******************************

At the MATLAB command prompt, type::

    >> matlabroot

    ans =

        '/Applications/MATLAB_R2019a.app'

.. warning::
    Ensure that your python environment has PhaRedox installed before continuing. To ensure that you do have it installed,
    execute the following code at your terminal: ``python -c "import pharedox ; print('success!')"``. You should see ``success!``
    at the terminal. If you see ``ModuleNotFoundError``, ensure that you (a) have installed pharedox and (b) are in the correct 
    virtual environment if you are using one.


In a system prompt, execute the following commands (replace ``<matlabroot>`` with
output of the above command)::

    $ cd <matlabroot>/extern/engines/python
    $ python setup.py install


How to set up a development environment
=======================================
If you'd like to work on the pipeline, you need to set up a few things on your computer
so you'll be as productive as possible.

Installing the development version of the code
**********************************************

Ensure that git is installed on your system. Git is software that helps us manage the
source code versions.

Once git is installed, download the source code. The easiest way is to issue the
following command::

    git clone https://github.com/ApfeldLab/pharynx_redox.git

This will download a new folder on your desktop called ``pharedox``.


Run the following command::

    pip install -e .

This will install PharRedox on your system such that when you change the source code,
the installed library itself will change. Finally, follow the instructions under the
`MATLAB`_ heading.

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

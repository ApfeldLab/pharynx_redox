.. _installation:

Installation
############

Basic Installation
==================

PhaRedox is a `Python <https://www.python.org/>`_ package. Many operating systems
ship with their own versions of Python, but these are either out of date or are
dangerous to mess with (as many system tools require them). Thus, it is recommended
that you install an alternate version of Python, separate from the one that ships
with your computer. I suggest that you use `Anaconda <https://www.anaconda.com/>`_ to
manage your Python versions. Anaconda is the most popular Python environment
management tool for scientists. However, if you prefer to manage your own Python
environments with `pyenv <https://github.com/pyenv/pyenv>`_ and `pyenv-virtualenv
<https://github.com/pyenv/pyenv-virtualenv>`_, that will work perfectly well too.
Anaconda is nice because it works the same in macOS, Windows, and Linux - which means
that I only have to write one set of installation instructions!

Installing Anaconda
-------------------

The easiest way to install Anaconda is to navigate to their `download page
<https://www.anaconda.com/products/individual>`_ then download and execute the
appropriate file for your operating system.

Creating a New Environment
--------------------------

The fastest way to create a new environment is through the terminal (macOS/Linux) or
the Anaconda Prompt (Windows, you can find this in your Start Menu). Simply execute
the following command::

    conda create --name pharedox python=3.7

This will create a new Python installation that is isolated from your operating
system's Python installation.

Installing the PhaRedox Library
-------------------------------

Now that you have a custom Python environment, we can install PhaRedox.

After activating your conda environment (`conda activate pharedox`), run::

    pip install pharedox

We're almost done!

MATLAB
~~~~~~

PhaRedox requires MATLAB for its 1D profile registration algorithm. Thus, we will
need to install that MATLAB library in your new Python environment. Unfortunately,
this process is difficult to automate - so it's left for you!

First, install MATLAB if you don't have it already.

.. warning::
    PhaRedox was developed using ``MATLAB_R2019a``. Other versions are not guaranteed
    to work, though you are free to try.

Open up MATLAB. Look for ``Set Path`` in the ``Home`` tab. Click it, then in the dialog,
click ``Add with subfolders``, and navigate to the PhaRedox source directory and select
the ``matlab`` folder. Finally, hit ``Save``.

Next, at the MATLAB command prompt, type::

    matlabroot

The system should print out a path that looks something like (on macOS)::

    '/Applications/MATLAB_R2019a.app'

On macOS/Linux, execute the following commands in Terminal. On Windows, execute them
in the Anaconda Prompt (you can find this in the Start Menu), and replace ``/`` with
``\``. In either case, replace ``<matlabroot>`` with the output of the above command.::

    cd <matlabroot>/extern/engines/python
    python setup.py install

Checking the Installation
-------------------------

If everything went well, you should have a working copy of PhaRedox, ready to analyze
your redox experiments! To make sure, activate your conda environment, and type the
following command::

    pharedox --help

You should see something like::

    Usage: pharedox [OPTIONS] COMMAND [ARGS]...

      Useful scripts for analyzing ratiometric microscopy data

    Options:
      --debug / --no-debug
      --help                Show this message and exit.

    Commands:
      analyze          Analyze an experiment
      create-settings  Create a settings file using the default template and...
      split-nc         Split an xarray DataArray into multiple Tiffs


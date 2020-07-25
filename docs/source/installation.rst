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

Installing Git
++++++++++++++

For the time being, we will install PhaRedox through its source code. The first step
is to download the source code via Git (the software that we use to track code
changes). By downloading PhaRedox through Git, you will be able to easily pull down
the latest changes.

In macOs, the easiest way to install Git is to try to run the following command::

    git --version

If you have Git installed, you're good to go. Otherwise, the system will prompt you
to install Xcode Command Line Tools. Simply follow the prompts and try the command
again - this time it should succeed.

In Windows, the easiest way is to install `GitHub Desktop <https://desktop.github
.com/>`. Simply download and run the installer.

Downloading the Source Code
+++++++++++++++++++++++++++

Now that Git is installed, we can download the source code. In macOS / Linux, open
your terminal, and navigate to where you would like the source code folder to reside::

    cd Path/to/parent/directory

Then, run the following command to download the source code into a new folder::

    git clone https://github.com/ApfeldLab/pharynx_redox.git

In Windows, follow `these <https://docs.github
.com/en/desktop/contributing-and-collaborating-using-github-desktop/cloning-a
-repository-from-github-to-github-desktop>`_ instructions. The link to the repository
(where you will find the green button in the instructions) is here: https://github
.com/ApfeldLab/pharynx_redox.

Installing the Code
+++++++++++++++++++

Activate your conda environment. In macOS, execute these commands in Terminal. In
Windows, execute them in Anaconda Prompt (you can find this in the Start Menu)::

    conda activate pharedox

Then, install the PhaRedox library::

    pip install -e Path/to/parent/directory/pharedox

In Windows, the slashes go the opposite way::

    pip install -e C:Path\to\parent\directory\pharedox

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



Updating PhaRedox
=================

To update PhaRedox, all you need to do is pull the latest changes from GitHub. On
macOS, execute the following commands::

    cd path/to/pharedox/
    git pull

On Windows, open GitHub Desktop, and follow the instructions `here <https://docs
.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/syncing
-your-branch#update-your-local-branch>`_ under ``Update your local branch``.

.. warning::
    If you made local changes to the code, this command might fail. This can happen if
    the changes that you made locally are different from the changes on the server. If
    this happens, we can force the remote changes to overwrite your local changes.
    Open a command prompt (on Windows, this can be done through the GitHub Desktop
    interface), and execute the following commands: ``git reset --hard HEAD`` and then
    ``git pull``. This should overwrite any of your local changes. Be careful!
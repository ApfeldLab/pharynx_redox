.. _installation:

############
Installation
############

This package requires at least Python 3.7. To check if your version is correct, execute
the following command::

    $ python --version
    Python 3.7.3

If the result is instead `Python 2.x.x`, or `Python 3.6.x` etc. you will need to
install a more up-to-date version of Python. `See here  <https://docs.python-guide
.org/starting/installation/>`_ for instructions on how to best do this.

Once you verify that the correct version of Python is installed, installing the
actual package is simple. First, you will need to download the source code. If git is
installed on your system this can be accomplished as follows::

    $ git clone https://github.com/half-adder/pharynx_redox.git

Otherwise, download and extract the following zip file: `<https://github
.com/half-adder/pharynx_redox/archive/master.zip>`_.

Once you have the folder, cd into the folder and install via pip (pip comes with
Python 3 by default)::

    $ cd pharynx_redox/
    $ python setup.py install


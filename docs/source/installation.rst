.. _installation:

############
Installation
############

Ensure that you have python installed on your system.::

    $ python --version
    Python 3.7.6

.. warning::
    Pharedox is not tested for Python < 3.7

Before moving forward, I would suggest that you use a virtual environments. Virtual environments keep packages that 
pharedox needs separate from your system's global python environment. See the instructions `here <https://python-guide-cn.readthedocs.io/en/latest/dev/virtualenvs.html>`_.
Note that this is *not required*, but is helpful in the long run.

If you think that you will want to edit the source code, follow the instructions on `how to set up your development environment`_. Otherwise,
after you have set up a virtual environment (or if you skipped that), simply execute the following line::

    pip install pharedox

MATLAB
======

If you would like to use the `1D Profile Registration` features, you will need to set up the MATLAB python engine. Ensure that MATLAB is
installed before continuing (your university probably has instructions on how to do this).

.. warning::
    PhaRedox was developed using ``MATLAB_R2019a``. Other versions are not guaranteed
    to work, though you are free to try.


Add PhaRedox files to MATLAB path
*********************************
In MATLAB, look for ``Set Path`` in the ``Home`` tab. Click it, then in the dialog,
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
    $ conda activate pharedox
    $ python setup.py install
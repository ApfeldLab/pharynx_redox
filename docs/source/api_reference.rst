#############
API reference
#############

This page provides an auto-generated summary of the APIs used for pharynx analysis.

Experiment
==========

.. py:currentmodule:: pharynx_redox.experiment

.. automodule:: pharynx_redox.experiment

.. autoclass:: Experiment

.. autoclass:: PairExperiment

.. autoclass:: CataExperiment


Pharynx IO
==========
This module coordinates loading data from and saving data to disk.

.. py:currentmodule:: pharynx_redox.pharynx_io

.. autosummary::
    :toctree: autosummary

    load_tiff_from_disk
    save_images_xarray_to_disk
    process_imaging_scheme_str
    load_images
    save_split_images_to_disk
    load_strain_map_from_disk
    load_all_rot_fl
    load_all_rot_seg

.. automodule:: pharynx_redox.pharynx_io
    :members:

Image Processing
================
This module contains the code for processing images, including segmentation, midline
calculation and measurement, translation, binary morphological operations, etc.

.. automodule:: pharynx_redox.image_processing
    :members:

Profile Processing
==================
This module contains the code for processing the measured profiles, including functions
for transforming ratios to OxD and E, trimming profiles, and registering profiles.

.. automodule:: pharynx_redox.profile_processing
    :members:

Plotting
========
This module contains plotting code

.. automodule:: pharynx_redox.plots
    :members:

Utils
=====

.. automodule:: pharynx_redox.utils
    :members:

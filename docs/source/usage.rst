###########
Basic Usage
###########

This package makes analysis of your image stacks easy. Here, we will walk through a
basic example of analyzing a fresh image stack from start to finish.

Before beginning, ensure that the package is successfully installed (see
:ref:`installation`).

Capturing Images
================

During image capture, be sure to keep track of the strain of each animal. You will
need to know who is who during analysis. This is referred to as the ``strain indexer``.
The software expects the strain indexer in the following format:

    +--------+--------------+------------+
    | strain | start_animal | end_animal |
    +--------+--------------+------------+
    | HD233  | 1            | 100        |
    +--------+--------------+------------+
    | SAY47  | 101          | 200        |
    +--------+--------------+------------+
    | HD233  | 201          | 300        |
    +--------+--------------+------------+

Organizing Files
================

Create a folder on your computer with the following format (we will refer to this as
the ``experiment ID``)::

    YYYY_MM_DD-experiment_identifiers

For example, let's say we've imaged HD233 and SAY47 at 4mm levamisole. An appropriate
experiment ID might be::

    2019_02_26-HD233_SAY47_4mm_lev

Once the folder is created, save your image stack as ``experiment_id.tiff`` and
strain indexer as ``experiment_id-indexer.csv`` in that folder. So in our example, we
would save our image stack as ``2019_02_26-HD233_SAY47_4mm_lev.tiff`` and our strain
indexer as ``2019_02_26-HD233_SAY47_4mm_lev-indexer.csv``.

Running the Analysis
====================

Once all of the files are in place, running the analysis is easy.

TODO: write run documentation
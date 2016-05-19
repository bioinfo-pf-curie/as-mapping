.. _QS:

.. Kenzo-Hugo Hillion
.. PROJECT_NAME
.. VERSION
.. 2016-05-18


Quick Start Guide for single-end reads Simulation
*****************

This page is a quick start guide, please read the full `online manual <link>`_ for more information.

See NEWS for information about changes in this and previous versions

What is it ?
============

| This pipeline was designed for the simulation of single-end reads

For any question about the pipeline, please contact kenzo.hillion@curie.fr

How to install it ?
===================

The following dependancies are required :

* BEDTools (version 2.21.0)
* vcf2diploid (version 0.2.6a)
* samtools (version )
* `ART <http://www.niehs.nih.gov/research/resources/software/biostatistics/art/>`_ (version ART-ChocolateCherryCake-03-19-2015) 

| To install :

.. code-block:: guess

  tar -zxvf master.tar.gz
  cd pip-master
  make CONFIG_SYS=config-install.txt install


+---------------+-----------------------------------------------------------------------------+
| SYSTEM CONFIGURATION                                                                        |
+===============+=============================================================================+
| PREFIX        | Path to installation folder                                                 |
+---------------+-----------------------------------------------------------------------------+
| ...           | ...                                                                         |
+---------------+-----------------------------------------------------------------------------+


Input Files
================

1. **A BED file** of ...

::

   chr1   0       16007
   chr1   16007   24571
   (...)

(...)


How to use it ?
===============

1. First have a look at the help message !

.. code-block:: guess

  command --help
  usage : command -h
  Use option -h|--help for more information

(...)

Test dataset
============


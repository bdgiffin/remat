.. REMAT documentation master file, created by
   sphinx-quickstart on Sun Nov 12 13:23:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/remat_logo.png
   :width: 640

*RE*\versible Physics and *MAT*\erial Model Library
===================================================

**Author**: `Brian Doran Giffin <https://github.com/bdgiffin>`_ (`brian.giffin@okstate.edu <mailto:brian.giffin@okstate.edu>`_)

REMAT is a proof-of-concept reversible physics library predicated on the use of fixed-precision numerics and integer arithmetic operations to ensure that all computations are exactly bit-reversible. Consequently, prior solution states may be precisely *rematerialized* through direct reversal of the forward-in-time operations. This is useful for the purpose of computing adjoint sensitivities of a time-dependent simulation without the need for solution checkpointing.

The core functionality of REMAT is primarily written in C++, with a Python API wrapper module provided to facilitate the use of REMAT in Python projects. Several interactive :ref:`Examples` are provided to demonstrate some of the core capabilities of REMAT.

Information regarding the :ref:`Installation` of REMAT and additional :ref:`Background` may be found within the accompanying documentation pages.

.. note::

   The documentation for REMAT is under active development. Please check back periodically for updates.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home <self>
   installation
   organization
   examples/index
   background

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

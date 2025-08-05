.. _Organization:

Organization
============

REMAT is primarily written as a C++ library, but comes equipped with a Python API to facilitate integration with other packages and application-driven workflows. The provided :ref:`Examples` offer an illustration of how these modules can be used to create and run various problems.

The accompanying Python modules are separated by their provided functionality, as follows:

 - ``REMAT.py``: The main API module providing Python bindings and wrapper utilities to interface with REMAT's underlying C++ framework.
 - ``Animation.py``: A graphical UI module that enables real-time interactive simulations using REMAT, primarily intended for demonstration and debugging purposes.
 - ``ExodusIO.py``: An output module enabling the export of REMAT simulation data to the Exodus file format.

Additionally, several Python-based pre-processing modules are included with REMAT to facilitate mesh generation and problem setup:

 - ``GeometryFactory.py``: A mesh management module that enables simple meshing, and importing mesh data from `Exodus <https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#exodus>`_, `DXF <https://help.autodesk.com/view/OARX/2024/ENU/?guid=GUID-235B22E0-A567-4CF6-92D3-38A2306D73F3>`_ (for truss/frame structures), and `meshio <https://github.com/nschloe/meshio>`_ sources (compatible with mesh data produced by `pygmsh <https://pypi.org/project/pygmsh/>`_).
 - ``Model.py``: A model management module, allowing for the definition of multiple mesh Parts, boundary and initial conditions, and contact interactions.

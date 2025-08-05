.. _Installation:

Installation
============

To build REMAT, you will first need to clone the `remat repository <https://github.com/bdgiffin/remat>`_ and build the underlying C++ shared object library files from source. If you are cloning a new repository, you will need to obtain all required submodules (BLT) via:

   .. code-block::

      git submodule update --init --recursive

To configure, build, and locally install REMAT:

   .. code-block::

      mkdir build
      cd build
      cmake ..
      make install

This will create a new local ``install`` subdirectory within the root REMAT project directory, containing all packaged Python modules and compiled C++ shared object libraries required to import and run the Python-based :ref:`Examples`.

If compiling the REMAT library to WebAssembly, be sure to configure CMake to use the Emscripten toolchain using ``emcmake`` to wrap the call to ``cmake``:

   .. code-block::

      emcmake cmake ..

Otherwise, the build/installation process described above remains the same.

To import the locally installed REMAT package within your Python project:

   .. code-block::

      import sys
      sys.path.append("PATH/TO/REMAT/install/package/")
      import REMAT

The ``REMAT/examples`` subdirectory provides further illustrative cases of invocations of the library within the context of a Python workflow.

When packaging a pygame project that uses REMAT using pygbag for execution in a web browser, make sure that you have compiled the REMAT library using Emscripten (see above note on configuring CMake using ``emcmake``), and include the ``install/package`` files directly within the local directory for your project. The provided ``examples`` offer a demonstration of how this can be accomplished, with the appropriate invocations of pygbag included in the ``examples/Makefile``.

Dependencies
------------

The core functionality of REMAT is primarily written in C++. A Python API wrapper module (implemented using `ctypes <https://docs.python.org/3/library/ctypes.html>`_) is provided to facilitate the use of REMAT in Python projects. Additional pre- and post-processing utilities are provided in Python, and demonstrated through several accompanying examples. While the underlying C++ framework does not have any dependencies apart from the standard library, the following Python packages are used by the API wrapper module and pre/post-processing utilities:

 - `numpy <https://numpy.org>`_
 - `scipy <https://scipy.org>`_
 - `pygmsh <https://pypi.org/project/pygmsh/>`_
 - `ezdxf <https://ezdxf.mozman.at/docs/>`_
 - `pyexodus <https://salvushub.github.io/pyexodus/>`_

REMAT utilizes a `CMake <https://cmake.org>`_-based build framework in combination with the funtionality provided by `BLT <https://computing.llnl.gov/projects/blt-build-link-test>`_. CMake version 3.14+ is required to build/install REMAT from source, whereas BLT is included as a git submodule.

The accompanying :ref:`Examples` additionally make use of the following Python packages to enable interactive visualization:

 - `pygame <https://pypi.org/project/pygame/>`_
 - `pygame_widgets <https://pygamewidgets.readthedocs.io/en/latest/>`_

Many of the :ref:`Examples` can be built for execution in a web browser using `Emscripten <https://emscripten.org>`_ to compile the underlying C++ library into WebAssembly, and `pygbag <https://pypi.org/project/pygbag/>`_ to package the pygame-based visualizations.

# REMAT

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14851878.svg)](https://doi.org/10.5281/zenodo.14851878)

**Documentation:** https://bdgiffin.github.io/remat/

---

REMAT is a proof-of-concept reversible physics library predicated on the use of fixed-precision numerics and integer arithmetic operations to ensure that all computations are exactly bit-reversible. Consequently, prior solution states may be precisely *rematerialized* through direct reversal of the forward-in-time operations. This is useful for the purpose of computing adjoint sensitivities of a time-dependent simulation without the need for solution checkpointing.

---

A central concept underpinning the developed bit-reversible computations pertains to the use of "dual" numbers which preserve information (bits) that would otherwise be lost due to finite precision rounding errors. In particular, consider the following example involving integer addition:

$$x=7 \qquad y=x+5=12$$

This operation may be exactly reversed through integer subtraction to recover the addend:

$$x=y-5=7$$

However, the operation of integer division cannot be exactly reversed through multiplication, i.e.

$$x=7 \qquad y=x\div5=1$$

and

$$x\neq y\times5=5$$

To ensure that the integer division operation can be exactly reversed, we must introduce an auxiliary variable to store the remainder:

$$x=7 \qquad y=x\div5=1 \qquad r=x \\, \text{mod} \\, 5=2$$

and thus the division may be reversed through the following correction involving the remainder:

$$x=y\times5+r=7$$

Ordinarily, the remainders of all such multiplication operations would need to be stored until the operation is later reversed, but this may entail intractable memory requirements if many such operations are carried out. To reduce memory overheads, a similar approach to that of [Maclaurin (2015)](https://proceedings.mlr.press/v37/maclaurin15.html) is adopted, whereby the dual variables $x^{\*}$ and $y^{\*}$ are introduced such that:

$$x=7 \qquad x^{\*}=n$$

$$y=x\div5=1 \qquad y^{\*}=x^{\*}\times5+x \\, \text{mod} \\, 5=n\times5+2$$

and

$$x=y\times5+y^{\*} \\, \text{mod} \\, 5=7 \qquad x^{\*}=y^{\*}\div5=n$$

This concept is applied such that any persistent stable variables $x$ are endowed with corresponding dual (adjoint) state variables $x^{\*}$ subject to an extended set of rules for all basic arithmetic operations which ensure that the pairing of $x$ and $x^{\*}$ efficiently preserves all bits that would otherwise be erased due to round-off.

---

## Dependencies

The core functionality of REMAT is primarily written in C++. A Python API wrapper module (implemented using [ctypes](https://docs.python.org/3/library/ctypes.html)) is provided to facilitate the use of REMAT in Python projects. Additional pre- and post-processing utilities are provided in Python, and demonstrated through several accompanying examples. While the underlying C++ framework does not have any dependencies apart from the standard library, the following Python packages are used by the API wrapper module and pre/post-processing utilities:
 - [numpy](https://numpy.org)
 - [scipy](https://scipy.org)
 - [pygmsh](https://pypi.org/project/pygmsh/)
 - [ezdxf](https://ezdxf.mozman.at/docs/)
 - [pyexodus](https://salvushub.github.io/pyexodus/)

REMAT utilizes a [CMake](https://cmake.org)-based build framework in combination with the funtionality provided by [BLT](https://computing.llnl.gov/projects/blt-build-link-test).

The accompanying `examples` additionally make use of the following Python packages to enable interactive visualization:
 - [pygame](https://pypi.org/project/pygame/)
 - [pygame_widgets](https://pygamewidgets.readthedocs.io/en/latest/)

Many of the `examples` can be compiled for execution in a web browser using [Emscripten](https://emscripten.org) to compile the underlying C++ library into WebAssembly, and [pygbag](https://pypi.org/project/pygbag/) to package the pygame-based visualizations.

---

## Getting started

To build REMAT, you will first need to clone this repository and build the underlying C++ shared object library files from source. If you are cloning a new repository, you will need to obtain all required submodules via:
```
git submodule update --init --recursive
```

To configure, build, and locally install REMAT:
```
mkdir build
cd build
cmake ..
make install
```
This will create a new local `install` subdirectory within the root REMAT project directory, containing all packaged Python modules and compiled C++ shared object libraries required to import and run the Python-based `examples`.

If compiling the REMAT library to WebAssembly, be sure to configure CMake to use the Emscripten toolchain using `emcmake` to wrap the call to `cmake`:
```
emcmake cmake ..
```
Otherwise, the build/installation process described above remains the same.

To import the locally installed REMAT package within your Python project:
```
import sys
sys.path.append("PATH/TO/REMAT/install/package/")
import REMAT
```
The `examples` subdirectory provides further illustrative cases of invocations of the library within the context of a Python workflow.

When packaging a pygame project that uses REMAT using pygbag for execution in a web browser, make sure that you have compiled the REMAT library using Emscripten (see above note on configuring CMake using `emcmake`), and include the `install/package` files within the local directory for your project. The provided `examples` offer a demonstration of how this can be accomplished, with the appropriate invocations of pygbag included in the `examples/Makefile`.
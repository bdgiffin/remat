# ReMat

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14851878.svg)](https://doi.org/10.5281/zenodo.14851878)

ReMat is a proof-of-concept reversible physics library predicated on the use of fixed-precision numerics and integer arithmetic operations to ensure that all computations are exactly bit-reversible. Consequently, prior solution states may be precisely *rematerialized* through direct reversal of the forward-in-time operations. This is useful for the purpose of computing adjoint sensitivities of a time-dependent simulation without the need for solution checkpointing.

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

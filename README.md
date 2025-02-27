# ECDSA Playground

## About

A collection of ECDSA implementations, typically found on resource-constrained devices such as smartcards.
These implementations are intended to analyze, test and play with them, e.g. when reversing
black box implementations.

These implementations are _not_ ready for production use! In particular, they are not side-channel resistant, although some
feature implementation details that are intended to study side-channel resistance.

## Implementations

- `01_plain_double_add`: textbook ECC *double & add* implemented with *affine* coordinates. Plain ECDSA.
  
- `02_projective_double_add`: textbook *double & add* implemented with *projective* coordinates. Scalar blinding
  with user-provided random. Scalar is not blinded.
    
- `03_projective_comb`: *double & add always* implementation with *projective* coordinates. 
  Scalar blinding with user-provided random. Scalar multiplication with comb-method
  and window-size 2. The scalar itself is not blinded. For further 
  reference cf. [1]
  
- `04_shamir_trick` : *double & add* using *projective* coordinates and *Shamir's trick*. Shamir's trick
  is here not used for signature verification, but rather for scalar multiplication by splitting
  kP = k_0P + k_1 \lambda P, where k_0 and k_1 are the lower resp. upper half of bits of k. \lambda P
  is precomputed. For further reference cf. [2] and [3]
  - variant a): the scalar is represented as a binary number. Since the case (0,0) does not require
    a point addition, this is per-se not sca-resistant (This is why `03_projective_comb` adds dummy operations).
  - variant b): the scalar is represented in ternary form with minimal hamming weight (non-adjacent form)
    Having several zeros improves performance, since fewer point additions are required. On the other hand,
    this again is per-se not sca-resistant, since we have different timings for zero/non-zero values.

- `05_montgomery` : straight-forward montgomery implementation using montgomery multiplication; nor further
   side-channel resistance

- `05_babygiantstep` : A toy implementation to illustrate the baby-step giant step algorithm. Of course, 
   this does not scale in python for nistp256, thus here we limit ourselves to 32 bit scalars.

## PyTest

Some tests are implemented to make sure that the implementations generate valid ECDSA signatures. To run
these tests, install `pytest`, open a commandline, change to the appropriate subdirectory and run `python -m pytest`, e.g. 

````
cd ECC_01_plain_double_add
python -m pytest
````
  
## Paper References

- [1] Victor Lomne and Thomas Roche: A Side Journey to Titan
- [2] Jerome A.: Low-Weight Binary Representations for Pairs of Integers
- [3] Franck Rondepierre: Revisiting Atomic Patterns for Scalar Multiplications on Elliptic Curves


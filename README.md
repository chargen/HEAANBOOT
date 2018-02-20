# HEAANBOOT
HEAANBOOT is software library that implements homomorphic encryption (HE) that supports fixed point arithmetics.
This library supports approximate operations between rational numbers.
The approximate error depends on some parameters and almost same with floating point operation errors.
The scheme in this library is on the paper "Homomorphic Encryption for Arithmetic of Approximate Numbers" (https://eprint.iacr.org/2016/421.pdf).

This library contain bootstrapping.

The bootsrapping technique is on the paper "Bootstrapping for Approximate Homomorphic Encryption" (https://eprint.iacr.org/2018/153.pdf).

This library is written by c++ and using NTL library (http://www.shoup.net/ntl/).

## How to use this library?
### 1. Build a static library/Running test functions
You can make a static library by typing make all in the /lib directory. After successful compilation you can find a static library libHEAANBOOT.a in the /lib directory.

You can run a test program. You need uncomment tests you need in /src/HEAAN.cpp and type make all in the /run directory. After successful compilation, you can see HEAAN file in the /run directory. To run the code type ./HEAANBOOT

We checked the program was working well on Ubuntu 16.04.2 LTS. You need to install NTL (with GMP), pThread, libraries. 

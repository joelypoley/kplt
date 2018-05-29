# KPLT

This repository contains my personal notes and some code I wrote while reading *[On the quaternion ell-isogeny problem](https://arxiv.org/abs/1406.0981)* by David Kohel, Kristin Lauter, Christophe Petit, Jean-Pierre Tignol.

# If you are going to clone this repo
Install `git lfs` first.
Instructions are [here](https://github.com/git-lfs/git-lfs/wiki/Installation).

`git lfs` is used to track the pdfs.

# Run tests
```bash
$ cd kplt
$ sage -python src/kplt_test.py
__main__.IdealsTest.test_connecting_ideal: 0.323
__main__.IdealsTest.test_element_of_norm: 0.007
__main__.IdealsTest.test_element_of_norm_large: 0.004
__main__.IdealsTest.test_prime_norm_representative: 0.003
__main__.IdealsTest.test_special_ell_power_equiv: 0.076
__main__.IdealsTest.test_strong_approximation: 0.037
----------------------------------------------------------------------
Ran 6 tests in 0.450s

OK
$ sage -python src/kplt_stress_test.py
p =  1000000007 I =  Fractional ideal (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)
p =  1000000007 I =  Fractional ideal (1/2 + 278800001943/2*j + 104550000730*k, 1/2*i + 69700000491*j + 278800001943/2*k, 174250001221*j, 174250001221*k)
p =  1000000007 I =  Fractional ideal (1/2 + 13888889/2*j + 6944445*k, 1/2*i + 6944444*j + 13888889/2*k, 13888889*j, 13888889*k)
...
```

# Files

## kplt.tex
This is the main TeX file, it describes the algorithm.

## quat_computing.tex
This file contains information on how you compute with quaternion algebras.

## src/
This folder contains the code.

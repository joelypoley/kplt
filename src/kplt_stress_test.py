from __future__ import print_function

import time

from sage.all import *

from kplt import prime_norm_representative
from kplt import element_of_norm
from kplt import left_ideal
from kplt import strong_approximation
from kplt import special_ell_power_equiv
from kplt import solve_ideal_equation
from kplt import connecting_ideal
from kplt import ell_power_equiv

large_primes = (next_prime(x) for x in range(1000000000, 1010000000))
primes_generator = (x for x in large_primes if mod(x, 4) == 3)
for p in primes_generator:
    B = QuaternionAlgebra(p)
    O = B.maximal_order()
    ell = Integer(2)
    alpha = O.random_element()
    if alpha.reduced_norm() == 0:
        continue
    I = left_ideal(
        [alpha, choice(Integer(alpha.reduced_norm()).divisors())], O)
    if I == O: print('I == O')
    print('p = ', p, 'I = ', I)
    if p == 311: continue
    start = time.time()
    # The assert statements at in the function ell_power_equiv ensure that the 
    # result is correct.
    _ = ell_power_equiv(I, O, ell)
    end = time.time() - start
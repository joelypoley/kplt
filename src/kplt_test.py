from __future__ import print_function

import unittest

from sage.all import *

from kplt import prime_norm_representative
from kplt import element_of_norm
from kplt import find_generators
from kplt import left_ideal
from kplt import strong_approximation 
from kplt import ell_power_equiv 
from kplt import solve_ideal_equation
from kplt import connecting_ideal

set_random_seed(0)

class IdealsTest(unittest.TestCase):

    def test_connecting_ideal(self):
        B = QuaternionAlgebra(59)
        i, j, k = B.gens()
        O_1 = B.maximal_order()
        O_2 = B.quaternion_order([(1 + k) / 2, (i + j) / 2, j, k], check=True)
        J = connecting_ideal(O_1, O_2)
        self.assertTrue(J.left_order() == O_1)
        self.assertTrue(J.right_order() == O_2)
        self.assertTrue(all(x in O_1 for x in J.basis()))

    def test_prime_norm_representative(self):
        B = QuaternionAlgebra(59)
        O = B.maximal_order()
        I = O.left_ideal(O.basis()).scale(2)
        J = prime_norm_representative(I, O, 4, 2)
        self.assertTrue(J.norm() in Primes() and J.left_order() == O)

    def test_find_generator(self):
        B = QuaternionAlgebra(59)
        O = B.maximal_order()
        I = O.left_ideal([2*x for x in O.basis()])
        N, alpha = find_generators(I)
        J = left_ideal([N, alpha], O)
        self.assertTrue(J == I)

    def test_strong_approximation(self):
        B = QuaternionAlgebra(59)
        ell = 3
        O = B.maximal_order()
        i, j, k = B.gens()
        gamma = 16 + 69*i
        alpha = 10*j + 2*k
        N = 101
        I = left_ideal([N, alpha], O)
        D = 4
        mu_0 = solve_ideal_equation(gamma, I, D, N, O)
        mu = strong_approximation(mu_0, N, O, ell)
        self.assertTrue(Integer(mu.reduced_norm()).prime_factors() == [ell])

    def test_element_of_norm(self):
        B = QuaternionAlgebra(59)
        O = B.maximal_order()
        I = O.left_ideal(O.basis()).scale(2)
        M = 200001
        gamma = element_of_norm(M, O)
        self.assertTrue(gamma.reduced_norm() == M)

    def test_element_of_norm_large(self):
        B = QuaternionAlgebra(1019)
        O = B.maximal_order()
        I = O.left_ideal(O.basis()).scale(2)
        M = 200000
        gamma = element_of_norm(M, O)
        self.assertTrue(gamma.reduced_norm() == M)

    def test_ell_power_equiv(self):
        B = QuaternionAlgebra(59)
        ell = Integer(3)
        O = B.maximal_order()
        i, j, k = B.gens()
        alpha = 1 + 2 * j
        I = left_ideal([alpha, 13], O)
        J = ell_power_equiv(I, O, ell)
        self.assertTrue(J.left_order() == O)
        self.assertTrue(Integer(J.norm()).prime_factors() == [ell])
        self.assertTrue([x in O for x in J.basis()])

if __name__ == '__main__':
    unittest.main()
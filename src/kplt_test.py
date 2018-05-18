from __future__ import print_function

import unittest

from sage.all import *
from kplt import prime_norm_representative
from kplt import element_of_norm
from kplt import find_generators
from kplt import left_ideal
from kplt import solve_congruence 
from kplt import strong_approximation 
from kplt import ell_power_equiv 

class IdealsTest(unittest.TestCase):

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

    def test_solve_congruence(self):
        B = QuaternionAlgebra(59)
        O = B.maximal_order()
        i, j, k = B.gens()
        gamma = 16 + 69*i
        alpha = 10*j + 2*k
        N = 101
        D = 4
        mu = solve_congruence(gamma, alpha, D, N, O)
        self.assertTrue(mu in B)
        NO = O.left_ideal(O.basis()).scale(N)
        self.assertTrue((gamma*mu - alpha) in NO)

    def test_strong_approximation(self):
        B = QuaternionAlgebra(59)
        ell = 3
        O = B.maximal_order()
        i, j, k = B.gens()
        gamma = 16 + 69*i
        alpha = 10*j + 2*k
        N = 101
        D = 4
        mu_0 = solve_congruence(gamma, alpha, D, N, O)
        mu, lamb = strong_approximation(mu_0, N, O, ell)
        self.assertTrue(Integer(mu.reduced_norm()).prime_factors() == [ell])
        NO = O.left_ideal(O.basis()).scale(N)
        self.assertTrue(mu - lamb*mu_0 in NO)

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
        alpha = Integer(1)/Integer(2) - (Integer(1)/Integer(2))*i + (Integer(3)/Integer(2))*j - (Integer(3)/Integer(2))*k
        N = 2300
        I = left_ideal([O.random_element(), 2300], O)
        J = ell_power_equiv(I, O, ell)
        self.assertTrue(J.left_order() == O)
        self.assertTrue(Integer(J.norm()).prime_factors() == [ell])
        self.assertTrue([x in O for x in J.basis()])

if __name__ == '__main__':
    unittest.main()
from __future__ import print_function

import unittest

from sage.all import *
from kplt import prime_norm_representative
from kplt import element_of_norm
from kplt import find_generators
from kplt import left_ideal
from kplt import solve_congruence 

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
        alpha = 10*j -k
        N = 101
        D = 4
        mu = solve_congruence(gamma, alpha, D, N, O)
        self.assertTrue(mu in B)
        NO = O.left_ideal(O.basis()).scale(N)
        self.assertTrue(Integer((gamma*mu - alpha) in NO))

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


if __name__ == '__main__':
    unittest.main()
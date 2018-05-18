from __future__ import print_function

import unittest

from sage.all import *
from kplt import prime_norm_representative
from kplt import element_of_norm


class IdealsTest(unittest.TestCase):

    def test_prime_norm_representative(self):
        B = QuaternionAlgebra(59)
        O = B.maximal_order()
        I = O.left_ideal(O.basis()).scale(2)
        J = prime_norm_representative(I, O, 4, 2)
        self.assertTrue(J.norm() in Primes() and J.left_order() == O)

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
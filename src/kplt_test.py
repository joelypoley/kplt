from __future__ import print_function

import unittest
import time

from sage.all import *

from kplt import prime_norm_representative
from kplt import element_of_norm
from kplt import left_ideal
from kplt import strong_approximation
from kplt import special_ell_power_equiv
from kplt import solve_ideal_equation
from kplt import connecting_ideal

set_random_seed(0)


class IdealsTest(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

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
        J, _ = prime_norm_representative(I, O, 4, 2)
        self.assertTrue(is_prime(Integer(J.norm())) and J.left_order() == O)

    def test_strong_approximation(self):
        B = QuaternionAlgebra(59)
        ell = 3
        O = B.maximal_order()
        i, j, k = B.gens()
        N = next_prime(100000)
        mu_0 = 16*j + 24 * k
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

    def test_special_ell_power_equiv(self):
        B = QuaternionAlgebra(59)
        ell = Integer(3)
        O = B.maximal_order()
        i, j, k = B.gens()
        alpha = 1 + 2 * j
        I = left_ideal([alpha, 13], O)
        J, _ = special_ell_power_equiv(I, O, ell)
        self.assertTrue(J.left_order() == O)
        self.assertTrue(Integer(J.norm()).prime_factors() == [ell])
        self.assertTrue([x in O for x in J.basis()])

    def ell_power_equiv(self):
        B = QuaternionAlgebra(59)
        ell = Integer(3)
        i, j, k = B.gens()
        gens = [(1 + k) / 2, (i + j) / 2, j, k]
        O = B.quaternion_order(gens)
        alpha = 2 * i - 2 * j + 2 * k
        assert alpha in O
        I = left_ideal([alpha, 24], O)
        J, _ = ell_power_equiv(I, O, ell)
        self.assertTrue(J.left_order() == O)
        self.assertTrue(Integer(J.norm()).prime_factors() == [ell])
        self.assertTrue([x in O for x in J.basis()])
        

if __name__ == '__main__':
    # We do this instead of unittest.main() so that the time for each test is 
    # printed.
    suite = unittest.TestLoader().loadTestsFromTestCase(IdealsTest)
    unittest.TextTestRunner(verbosity=0).run(suite)

from __future__ import print_function

import unittest
import time

from sage.all import *

from kplt import prime_norm_representative
from kplt import element_of_norm
from kplt import find_generators
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
        J = prime_norm_representative(I, O, 4, 2)
        self.assertTrue(J.norm() in Primes() and J.left_order() == O)

    def test_find_generator(self):
        B = QuaternionAlgebra(59)
        O = B.maximal_order()
        I = O.left_ideal([2 * x for x in O.basis()])
        N, alpha = find_generators(I)
        J = left_ideal([N, alpha], O)
        self.assertTrue(J == I)

    def test_strong_approximation(self):
        B = QuaternionAlgebra(59)
        ell = 3
        O = B.maximal_order()
        i, j, k = B.gens()
        gamma = 16 + 69 * i
        alpha = 10 * j + 2 * k
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

    def test_special_ell_power_equiv(self):
        B = QuaternionAlgebra(59)
        ell = Integer(3)
        O = B.maximal_order()
        i, j, k = B.gens()
        alpha = 1 + 2 * j
        I = left_ideal([alpha, 13], O)
        _, J = special_ell_power_equiv(I, O, ell)
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
        J = ell_power_equiv(I, O, ell)
        self.assertTrue(J.left_order() == O)
        self.assertTrue(Integer(J.norm()).prime_factors() == [ell])
        self.assertTrue([x in O for x in J.basis()])


class StressTest(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

    def test_prime_norm_representative(self):
        primes_generator = (x for x in Primes() if mod(x, 4) == 3)
        bound = 1000000
        for p in primes_generator:
            if p > bound: break
            B = QuaternionAlgebra(p)
            O = B.maximal_order()
            alpha = O.random_element()
            if alpha.reduced_norm() == 0: 
                continue
            I = left_ideal([alpha, choice(Integer(alpha.reduced_norm()).divisors())], O)
            J = prime_norm_representative(I, O, 4, 2)
            N = J.norm()
            if I == O: print('I == O')
            M = N * Integer(2)**50
            print('p = ', p, 'M = ', M)
            start = time.time()
            gamma = element_of_norm(M, O)
            end = time.time() - start 
            self.assertTrue(gamma.reduced_norm() == M)


    # def test_element_of_norm(self):
    #     B = QuaternionAlgebra(59)
    #     O = B.maximal_order()
    #     I = O.left_ideal(O.basis()).scale(2)
    #     M = 200001
        

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(IdealsTest)
    suite_2 = unittest.TestLoader().loadTestsFromTestCase(StressTest)
    unittest.TextTestRunner(verbosity=0).run(suite)
    unittest.TextTestRunner(verbosity=0).run(suite_2)

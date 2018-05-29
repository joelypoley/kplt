from __future__ import print_function

import unittest

from sage.all import *
from cornacchia import cornacchia


class CornacchiaTest(unittest.TestCase):

    def test_cornacchia(self):
        d = 2
        m = 22
        res = cornacchia(d, m)
        self.assertTrue(res is not None)
        x, y = res
        print(type(x))
        print(type(y))
        self.assertTrue(type(x) == type(y) == Integer)
        self.assertTrue(x ** 2 + d * y ** 2 == m)
        self.assertTrue(gcd(x, y) == 1)


if __name__ == "__main__":
    unittest.main()

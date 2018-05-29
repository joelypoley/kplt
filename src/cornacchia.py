from __future__ import print_function
from sage.all import *


def cornacchia(d, m):
    """Find a primitive solution to x^2 + d*y^2 = m.

    A primitive solution means that (x, y) = 1.

    Args:
        d: A sage integer.
        m: A sage integer.

    Returns:
        The soution (x, y) as a tuple of sage integers or None if there is no
        solution.
    """
    # print(m)
    if m <= 0:
        raise ValueError("m must be greater than 0 but m was " + str(m))

    if d == 0:
        raise ValueError("d must be nonzero")

    if not mod(-d, m).is_square():
        return None

    prev = m
    curr = Integer((mod(-d, m).sqrt()))
    while curr ** 2 >= m:
        prev, curr = curr, Integer(mod(prev, curr))

    x = curr
    y_squared = (m - x ** 2) / d

    if not y_squared in Integers():
        return None
    elif not Integer(y_squared).is_square():
        return None
    else:
        y = sqrt(Integer(y_squared))
        assert x ** 2 + d * y ** 2 == m and gcd(x, y) == 1
        return x, y

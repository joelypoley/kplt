from __future__ import print_function

from sage.all import *
from cornacchia import cornacchia


def normalized_norm(alpha, I):
    """
    Returns the normalized norm N(alpha) / N(I). Using the notation in the
    paper, it computes q_I(alpha).

    Args:
        alpha: An element of I.
        I: An fractional ideal.  

    Returns:
        An Integer.
    """
    return Integer(alpha.reduced_norm() / I.norm())


def one_norm(alpha):
    """Returns the 1-norm of alpha.
        Args:
            alpha: An element with a coefficient_tuple() method.

        Returns:
            ||alpha||_1.
    """
    return sum(abs(x_i) for x_i in alpha.coefficient_tuple())


def is_minkowski_basis(basis):
    """Checks if the basis is Minkowski."""
    # TODO: Ask Steven about Minkowskiness.
    return all(
        one_norm(alpha_i) <= one_norm(
            sum(alpha_j for alpha_j in basis if alpha_j != alpha_i))
        for alpha_i in basis)


def prime_norm_representative(I, O):
    """
    Given a maximal order `O` and a left `O`-ideal return another
    left `O` J in the same class, but with prime norm.

    Args:
        I: A left O-ideal.
        O: A maxmimal order in a quaternion algebra.
    Returns:
        Another left O-ideal in the same class with prime norm.

    """
    if not is_minkowski_basis(I.basis()):
        print('Warning: The ideal I does not have a minkowski basis'
              'precomputed and Sage can not do it for you.')
    m = 10  # TODO: Change this to a proper bound.
    B = I.quaternion_algebra()
    if mod(B.discriminant(), 4) != 3:
        raise NotImplementedError('The quaternion algebra must have'
                                  ' discriminant p = 3 mod 4')
    alpha = B(0)
    # Choose random elements in I until one is found with prime norm.
    while normalized_norm(alpha, I) not in Primes():
        alpha = sum(
            (ZZ.random_element(-m, m + 1) * alpha_i for alpha_i in I.basis()))

    # Now we have an element alpha of prime norm the ideal J = I*gamma has
    # prime norm where gamma = conjguate(alpha) / N(I).
    N = I.norm()
    gamma = alpha.conjugate() / N
    J_basis = [alpha_i * gamma for alpha_i in I.basis()]
    J = O.left_ideal(J_basis)
    return J


def f(x, y, q):
    """Computes the principal form f(x, y) in the paper.
    Args:
        x: A Sage integer.
        y: A Sage integer.  
        q: i^2 if p = 3 mod 4.

    Returns:
        A Sage integer.
    """
    return x**2 + q * y**2


def solve_norm_equation(q, r):
    """ Solves x^2 + q*y^2 = r.

    Args:
        r: A Sage integer.
        q: A Sage integer.

    Returns:
        Tuple of Sage integers (x, y) or None if there is no solution.
    """
    # At the moment this is very easy since I am assuming we q = 1. It gets
    # harder without this assumption.
    if q == 1:
        try:
            sol = two_squares(r)
            return sol
        except ValueError:
            return None
    else:
        raise NotImplementedError('q must be equal to 1 and was ' + str(q))


def element_of_norm(M, O):
    """Finds an element of B with norm M.

    Args:
        M: A sage integer.
        O: A maximal order in a quaternion algebra.

    Returns:
        gamma in B such that gamma.reduced_norm() == M
    """
    B = O.quaternion_algebra()
    if mod(B.discriminant(), 4) != 3:
        raise NotImplementedError('The quaternion algebra must have'
                                  ' discriminant p = 3 mod 4')
    i = B.gen(0)
    j = B.gen(1)
    q, p = [Integer(-i**2), Integer(-j**2)]
    m = 10000
    r = 0
    sol = None
    x_2 = ZZ.random_element(x=m)
    y_2 = ZZ.random_element(x=m)
    while True:
        x_2 = ZZ.random_element(-m, m + 1)
        y_2 = ZZ.random_element(-m, m + 1)
        r = M - p * f(x_2, y_2, q)

        if r not in Primes():  # Can replace with easily factorable.
            continue

        # The norm equation is N(x + iy) = x^2 + qy^2.
        sol = cornacchia(q, r)
        if sol is not None:
            break

    x_1, y_1 = sol
    alpha = x_1 + y_1 * i
    beta = x_2 + y_2 * i
    assert Integer((alpha + beta * j).reduced_norm()) == M
    return alpha + beta * j


from __future__ import print_function

from sage.all import *
from sage.algebras.quatalg.quaternion_algebra import quaternion_algebra_cython
from sage.rings.finite_rings.integer_mod import mod
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.arith.misc import is_prime
from sage.arith.misc import xgcd
from sage.arith.misc import gcd
from sage.arith.misc import two_squares
from sage.arith.functions import lcm
from sage.sets.primes import Primes
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from cornacchia import cornacchia
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.matrix.constructor import matrix


def ideal_product(I, J):
    """Returns IJ."""
    assert I.right_order() == J.left_order()
    # Don't really know why I have to do this.
    mat_1 = matrix([x.coefficient_tuple() for x in I.right_order().basis()])
    mat_2 = matrix([x.coefficient_tuple() for x in I.left_order().basis()])
    matcoeff = mat_2 * ~mat_1
    N = lcm(x.denominator() for x in matcoeff.coefficients())

    K = left_ideal([N * x * y for x in I.basis() for y in J.basis()],
                   I.left_order())

    assert K.left_order() == I.left_order()
    assert K.right_order() == J.right_order()
    return K


def connecting_ideal(O_1, O_2):
    """Returns an 0_1, O_2-connecting ideal.

    Args:
        O_1: A maximal order in a rational quaternion algebra.
        O_2: A maximal order in the same quaternion algebra.

    Returns:
        An ideal I that is a left O_1 ideal and a right O_2 ideal. Moreover I
        is a subset of O_1.
    """
    mat_1 = matrix([x.coefficient_tuple() for x in O_1.basis()])
    mat_2 = matrix([x.coefficient_tuple() for x in O_2.basis()])
    matcoeff = mat_2 * ~mat_1
    N = lcm(x.denominator() for x in matcoeff.coefficients())
    J = left_ideal([N * x for x in O_2.basis()] +
                   [N * x * y for x in O_1.basis() for y in O_2.basis()], O_1)

    assert J.left_order() == O_1
    assert J.right_order() == O_2
    assert all(x in O_1 for x in J.basis())

    return J


def left_ideal(gens, O):
    """Returns the left O-ideal generated by gens.

    If gens = [gens_0, ..., gens_n] then the ideal return is O*gens_0 + ... +
    O*gens_n.

    Args:
        gens: A list of elements of O.
        O: An order in a quaternion algebra.

    Returns:
        The left O-ideal generated by gens.
    """
    if not all(x in O for x in gens):
        raise ValueError('All the generators must be in O.')

    all_gens = [x * y for x in O.basis() for y in gens]
    B = O.quaternion_algebra()
    Z, d = (quaternion_algebra_cython.
            integral_matrix_and_denom_from_rational_quaternions(all_gens))
    H = Z.hermite_form(include_zero_rows=False)
    basis = (quaternion_algebra_cython.
             rational_quaternions_from_integral_matrix_and_denom(B, H, d))

    I = O.left_ideal(basis)

    assert all(x in O for x in I.basis())
    assert all (x * y in O for x in O.basis() for y in I.basis())
    assert I.left_order() == O
    return I


def random_combination(basis, bound=1000):
    """Return a random integer linear combination of the elements of basis."""
    return sum(ZZ.random_element(-bound, bound + 1) * x_i for x_i in basis)


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
    return all(
        one_norm(alpha_i) <= one_norm(
            sum(alpha_j for alpha_j in basis if alpha_j != alpha_i))
        for alpha_i in basis)


def prime_norm_representative(I, O, D, ell):
    """
    Given an order O and a left O-ideal I return another
    left O-ideal J in the same class, but with prime norm.

    This corresponds to Step 1 in the notes. So given an ideal I it returns
    an ideal in the same class but with reduced norm N where N != ell is
    a large prime coprime to both D and p, and ell is a quadratic
    nonresidue module N.


    Args:
        I: A left O-ideal.
        O: An order in a quaternion algebra.
        D: An integer.
        ell: A prime.
    Returns:
        Another left O-ideal in the same class with prime norm N. N will
        be coprime to both D and p, and ell will be a nonquadratic residue
        module N. 

    """
    if not is_minkowski_basis(I.basis()):
        print('Warning: The ideal I does not have a minkowski basis'
              ' precomputed and Sage can not do it for you.')

    nrd_I = I.norm()
    alpha = B(0)
    normalized_norm = Integer(alpha.reduced_norm() / nrd_I)
    # Choose random elements in I until one is found with norm N*nrd(I) where N
    # is prime.
    m_power = 3
    m = Integer(2)**m_power  # TODO: Change this to a proper bound.
    count = 0
    while not is_prime(normalized_norm) or normalized_norm.divides(
            D) or normalized_norm == ell or normalized_norm == p or mod(
                ell, normalized_norm).is_square():
        # Make a new random element.
        alpha = random_combination(I.basis(), bound=m)
        normalized_norm = Integer(alpha.reduced_norm() / nrd_I)


        # Increase the box we search in if we've been trying for too long. Note
        # this was just a random heuristic I came up with, it's not in the
        # paper.
        count += 1
        if count > 4 * m_power:
            m_power += 1
            m = Integer(2)**m_power
            count = 0

    # We now have an element alpha with norm N*nrd(I) where N is prime. The
    # ideal J = I*gamma has prime norm where gamma = conjugate(alpha) / nrd(I).
    J = I.scale(gamma)

    assert is_prime(Integer(J.norm()))
    assert not mod(ell, Integer(J.norm())).is_square()
    assert gcd(Integer(J.norm()), D) == 1
    return J


def solve_norm_equation(q, r):
    """ Solves x^2 + q*y^2 = r.

    Args:
        r: A Sage integer.
        q: A Sage integer.

    Returns:
        Tuple of Sage integers (x, y) or None if there is no solution.
    """
    if q == 1:
        try:
            sol = two_squares(r)
            assert sol[0]**2 + q* sol[1]**2 == r
            return sol
        except ValueError:
            return None
    else:
        sol = cornacchia(q, r)
        assert True if sol in None else sol[0]**2 + q* sol[1]**2 == r
        return sol


def element_of_norm(M, O, bound=100):
    """Finds an element of B with norm M.

    This corresponds to Step 3 of the algorithm in the notes.

    Args:
        M: A sage integer.
        O: A maximal order in a quaternion algebra.
        bound: The values 0 <= x_2, y_2 <= 100 will be tried. If no solution is
            found then the function returns None.

    Returns:
        gamma in B such that gamma.reduced_norm() == M or None if there is no
        solution in the box [0, bound]**2.
    """
    B = O.quaternion_algebra()
    if mod(B.discriminant(), 4) != 3:
        raise NotImplementedError('The quaternion algebra must have'
                                  ' discriminant p = 3 mod 4')
    a, b = B.invariants()
    i, j, _ = B.gens()
    q, p = -Integer(a), -Integer(b)
    for x_2 in range(bound+1):
        for y_2 in range(bound+1):
            r = M - p * (x_2**2 + q * y_2**2)

            if r not in Primes():  # Can replace with easily factorable.
                continue

            # The norm equation is N(x + iy) = x^2 + qy^2.
            sol = solve_norm_equation(q, r)
            if sol is not None:
                x_1, y_1 = sol
                alpha = x_1 + y_1 * i
                beta = x_2 + y_2 * i
                assert Integer((alpha + beta * j).reduced_norm()) == M
                return alpha + beta * j

    return None

    


def find_generators(I):
    """Finds N, alpha in I such that I = ON + O*alpha

    This corresponds to Step 2 of the algorithm in the notes.

    Args:
        I: An ideal.

    Returns:
        A tuple (N, alpha) where N is a Sage integer and alpha is in I.
    """
    O = I.left_order()
    B = I.quaternion_algebra()
    N = I.norm()
    alpha = I.basis()[0]
    while alpha.reduced_norm() == 0 or Integer(alpha.reduced_norm()).divides(
            N**2):
        alpha = random_combination(I.basis())

    J = left_ideal([N, alpha], O)
    assert J == I
    return N, alpha


def solve_ideal_equation(gamma, I, D, N, O):
    """Find mu_0 in Rj such that (O* gamma / NO)[mu_0] = I / NO.

    Args:
        gamma: An element of O.
        I: A left O-ideal.
        D: The index [O : R + Rj].
        N: The norm of I. Must be prime.
        O: A special order in a quaternion algebra.

    Returns:
        mu_0 in Rj such that 0 != gamma * mu_0 in I.
    """
    d, c, _ = xgcd(D, N)
    assert d == 1
    a, b = [Integer(x) for x in O.quaternion_algebra().invariants()]
    F = GF(N)
    # The suffix _ff means that this element is over a finite field, not the
    # rationals.
    B_ff = QuaternionAlgebra(F, a, b)
    i_ff, j_ff, k_ff = B_ff.gens()

    def phi(alpha):
        t, x, y, z = [
            Integer(D * c * alpha_i) for alpha_i in alpha.coefficient_tuple()
        ]
        return t + x * i_ff + y * j_ff + z * k_ff

    gamma_ff = phi(gamma)
    gamma_ff_mat = gamma_ff.matrix(action='left')

    I_basis_ff = [phi(alpha).coefficient_tuple() for alpha in I.basis()]
    lin_system = matrix(F, [gamma_ff_mat[2], gamma_ff_mat[3]] + I_basis_ff)
    sol = lin_system.left_kernel().basis()[0]
    y, z = sol[0], sol[1]
    mu_ff = y * j_ff + z * k_ff

    B = O.quaternion_algebra()
    i, j, k = B.gens()
    mu_0 = sum(
        Integer(coeff) * elem
        for coeff, elem in zip(mu_ff.coefficient_tuple(), [1, i, j, k]))

    assert 0 != mu_0
    assert gamma * mu_0 in I, "gamma = " + str(gamma) + " m_0 = " + str(mu_0) + " I = " + str(I) + " O = " + str(O)

    return mu_0


def strong_approximation(mu_0, N, O, ell):
    """Find mu in O with nrd(mu) = ell^e and mu = lambda * mu_0 mod NO

    Args:
        mu_0: An element of Rj.
        N: A prime.
        O: A special order.
        ell: A prime.

    Returns:
        mu in O with nrd(mu) = ell^e and mu = lambda * mu_0 mod NO.
    """
    ell = Integer(ell)
    N = Integer(N)
    #print(N)
    B = O.quaternion_algebra()
    p = B.discriminant()
    i, j, k = B.gens()
    t_0, x_0, y_0, z_0 = mu_0.coefficient_tuple()
    assert t_0 == x_0 == 0
    beta_0 = y_0 + z_0 * i
    e = 100 if (
        ~mod(p * Integer(beta_0.reduced_norm()), N)).is_square() else 101
    lamb = Integer(
        (ell**e * ~mod(p * Integer(beta_0.reduced_norm()), N)).sqrt())
    assert mod(ell**e, N) == mod(lamb**2 * Integer(mu_0.reduced_norm()), N)
    mu = None
    while mu is None:
        assert N.divides(ell**e - p * lamb**2 * Integer(beta_0.reduced_norm()))
        lhs = Integer(
            (ell**e - p * lamb**2 * Integer(beta_0.reduced_norm())) / N)
        z_1 = ZZ.random_element(0, N**2)
        b = Integer(lhs - p * lamb * 2 * z_0 * z_1)
        y_1 = b * ~mod(2 * Integer(y_0) * p * lamb, N**2)
        beta_1 = Integer(y_1) + Integer(z_1) * i

        assert mod(lhs, N**2) == mod(
            p * lamb * (beta_0 * beta_1.conjugate()).reduced_trace(), N**2)
        assert (N**2).divides(ell**e - p *
                              (lamb * beta_0 + N * beta_1).reduced_norm())
        r = (ell**e - p * (lamb * beta_0 + N * beta_1).reduced_norm()) / N**2
        assert r > 0
        # TODO: Ensure the r is postive and in the right range.

        sol = solve_norm_equation(1, r)
        if sol is not None:
            t_1, x_1 = sol
            alpha_1 = t_1 + x_1 * i
            mu_1 = alpha_1 + beta_1 * j
            mu = lamb * mu_0 + N * mu_1
            break

    assert mu - lamb * mu_0 in O.left_ideal(O.basis()).scale(N)
    assert mu.reduced_norm() == ell**e
    return mu


def special_ell_power_equiv(I, O, ell, print_progress=False):
    """Solve ell isogeny problem.

    Args:
        I: A left O-ideal.
        O: A special order in a quaternion algebra.
        ell: A prime.
        print_progress: True if you want to print progress.

    Returns:
        An ideal non-fractional ideal J in the same class as I that has ell
        power norm.
    """
    ell = Integer(ell)
    D = 4
    I = prime_norm_representative(I, O, D, ell)
    N, alpha = find_generators(I)
    gamma = element_of_norm(N * ell**50, O)
    mu_0 = solve_ideal_equation(gamma, I, D, N, O)
    mu = strong_approximation(mu_0, N, O, ell)
    beta = gamma * mu
    J = I.scale(beta.conjugate() / N)
    assert J.left_order() == O
    assert Integer(J.norm()).prime_factors() == [ell]
    assert [x in O for x in J.basis()]
    return beta, J


def ell_power_equiv(J, O, ell, print_progress=False):
    B = O.quaternion_algebra()
    assert is_prime(B.discriminant()) and mod(B.discriminant(), 4) == 3
    O_special = B.maximal_order()
    I = connecting_ideal(O_special, O)
    K = ideal_product(I, J)
    gamma_1, I_1 = special_ell_power_equiv(I, O_special, ell)
    gamma_2, I_2 = special_ell_power_equiv(K, O_special, ell)
    gamma = gamma_1.conjugate() * gamma_2 / Integer(I.norm())
    assert J.left_order() == O
    assert Integer(J.norm()).prime_factors() == [ell]
    assert [x in O for x in J.basis()]
    return J.scale(gamma.conjugate() / Integer(J.norm()))

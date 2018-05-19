from __future__ import print_function

from sage.all import *
from sage.algebras.quatalg.quaternion_algebra import quaternion_algebra_cython
from cornacchia import cornacchia


def left_ideal(gens, O):
    """Returns a left O ideal generated by gens as a left O-ideal."""
    if not all(x in O for x in gens):
        raise ValueError('All the generators must be in O')

    all_gens = [x * y for x in O.basis() for y in gens]
    B = O.quaternion_algebra()
    Z, d = (quaternion_algebra_cython.
            integral_matrix_and_denom_from_rational_quaternions(all_gens))
    H = Z.hermite_form(include_zero_rows=False)
    basis = (quaternion_algebra_cython.
             rational_quaternions_from_integral_matrix_and_denom(B, H, d))

    return O.left_ideal(basis)


def random_combination(basis, bound=1000):
    """Return a random integer linear combination of the elements of basis."""
    return sum(ZZ.random_element(-bound, bound + 1) * x_i for x_i in basis)


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


def prime_norm_representative(I, O, D, ell):
    """
    Given a maximal order O and a left O-ideal return another
    left O-ideal J in the same class, but with prime norm N.

    This corresponds to Step 1 in the notes. So given an ideal I it returns
    an ideal in the same class but with reduced norm N where N != ell is
    a large prime coprime to D and p and such that ell is a quadratic
    nonresidue module N.


    Args:
        I: A left O-ideal.
        O: A maximal order in a quaternion algebra.
        D: An integer.
        ell: A prime.
    Returns:
        Another left O-ideal in the same class with prime norm N. N will
        be coprime to D and ell and ell and ell will be a nonquadratic
        residue module N. 

    """
    # Check preconditions.
    if not is_minkowski_basis(I.basis()):
        print('Warning: The ideal I does not have a minkowski basis'
              'precomputed and Sage can not do it for you.')

    m = 1000  # TODO: Change this to a proper bound.
    B = I.quaternion_algebra()
    p = B.discriminant()
    if mod(p, 4) != 3 or not is_prime(p):
        raise NotImplementedError('The quaternion algebra must have'
                                  ' prime discriminant p == 3 mod 4')

    if O != B.maximal_order():
        raise NotImplementedError(
            'The order O must be a special order. I.e. O must satisfy'
            ' O.quaternion_algebra().maximal_order() == O.')

    N = I.norm()
    alpha = B(0)
    normalized_norm = Integer(alpha.reduced_norm() / N)
    # Choose random elements in I until one is found with prime norm.
    while not is_prime(normalized_norm) or normalized_norm.divides(
            D) or normalized_norm == ell or normalized_norm == p or mod(
                ell, normalized_norm).is_square():
        alpha = random_combination(I.basis(), bound=m)
        normalized_norm = Integer(alpha.reduced_norm() / N)

    # Now we have an element alpha of prime norm the ideal J = I*gamma has
    # prime norm where gamma = conjguate(alpha) / N(I).
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
        return cornacchia(q, r)


def element_of_norm(M, O):
    """Finds an element of B with norm M.

    This corresponds to Step 3 of the algorithm in the notes.

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
    m = 100  # TODO: Replace with proper bound.
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
        sol = solve_norm_equation(q, r)
        if sol is not None:
            break

    x_1, y_1 = sol
    alpha = x_1 + y_1 * i
    beta = x_2 + y_2 * i
    assert Integer((alpha + beta * j).reduced_norm()) == M
    return alpha + beta * j


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
    if not J == I:
        raise RuntimeError('Oh dear. The algorithm failed.')

    return N, alpha


def solve_congruence(gamma, alpha, D, N, O):
    """Solve gamma * mu = alpha mod NO for mu in Zj + Zk and mu != 0 mod N.

    This corresponds to Step 4 in the algorithm.

    Args:
        gamma: An element of the quaternion algebra.
        alpha: An element of the quaternion algebra.
        N: A prime number.
        O: A maximal order in the quaternion algebra containg Z<i, j>.

    Returns:
        A solution mu to the equation gamma * mu = alpha mod NO for mu in 
        Zj + Zk and mu != 0 mod N.
    """

    d, a, b = xgcd(D, N)
    assert d == 1
    a, b = [Integer(x) for x in O.quaternion_algebra().invariants()]
    F = GF(N)
    B_ff = QuaternionAlgebra(F, a, b)
    # The suffix _ff means that this element is over a finite field, not the
    # rationals.
    i_ff, j_ff, k_ff = B_ff.gens()

    def phi(alpha):
        t, x, y, z = [
            Integer(D * a * alpha_i) for alpha_i in alpha.coefficient_tuple()
        ]
        return t + x * i_ff + y * j_ff + z * k_ff

    gamma_ff = phi(gamma)
    alpha_ff = phi(alpha)

    alpha_ff_vec = matrix(F, alpha_ff.coefficient_tuple())
    gamma_ff_mat = gamma_ff.matrix(action='left')
    gamma_ff_mat = matrix(F, [gamma_ff_mat[2], gamma_ff_mat[3]])
    sol = gamma_ff_mat.solve_left(alpha_ff_vec)
    y, z = [sol[0][c] for c in [0, 1]]
    mu_ff = y * j_ff + z * k_ff
    assert gamma_ff * mu_ff == alpha_ff

    B = O.quaternion_algebra()
    i, j, k = B.gens()
    return sum(
        Integer(coeff) * elem
        for coeff, elem in zip(mu_ff.coefficient_tuple(), [1, i, j, k]))


def solve_ideal_equation(gamma, I, D, N, O):
    d, a, b = xgcd(D, N)
    assert d == 1
    a, b = [Integer(x) for x in O.quaternion_algebra().invariants()]
    F = GF(N)
    B_ff = QuaternionAlgebra(F, a, b)
    # The suffix _ff means that this element is over a finite field, not the
    # rationals.
    i_ff, j_ff, k_ff = B_ff.gens()

    def phi(alpha):
        t, x, y, z = [
            Integer(D * a * alpha_i) for alpha_i in alpha.coefficient_tuple()
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
    return sum(
        Integer(coeff) * elem
        for coeff, elem in zip(mu_ff.coefficient_tuple(), [1, i, j, k]))


def strong_approximation(mu_0, N, O, ell):
    """Find mu in O with nrd(mu) = ell^e and mu = lambda * mu_0 mod NO"""
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
        #z_1 = ZZ.random_element(0, N**2)
        z_1 = ZZ.random_element(0, N)
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
            assert mu.reduced_norm() == ell**e
            break

    assert mu - lamb*mu_0 in O.left_ideal(O.basis()).scale(N)
    return mu


def ell_power_equiv(I, O, ell, print_progress=False):
    """Solve ell isogeny problem"""
    ell = Integer(ell)
    D = 4
    I = prime_norm_representative(I, O, D, ell)
    N, alpha = find_generators(I)
    gamma = element_of_norm(N * ell**50, O)
    mu_0 = solve_ideal_equation(gamma, I, D, N, O)
    mu = strong_approximation(mu_0, N, O, ell)
    beta = gamma * mu
    J = I.scale(beta.conjugate() / N)
    return J

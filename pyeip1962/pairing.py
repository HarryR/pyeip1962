from .group import AbstractPointG2, AbstractPointG1, AbstractGroup
from math import log2, floor


def linefunc(P1, P2, T):
    """
    Create a function representing the line between P1 and P2, and evaluate it at T
    tangent if P1 == P2.
    vp3 = vertical line through P3 = P1 + P2
    """
    assert P1 and P2 and T  # No points-at-infinity allowed, sorry
    x1, y1 = P1
    x2, y2 = P2
    xt, yt = T
    if x1 != x2:
        # Addition
        m = (y2 - y1) / (x2 - x1)
        return m * (xt - x1) - (yt - y1)
    elif y1 == y2:
        # Doubling
        m = 3 * x1**2 / (2 * y1)
        return m * (xt - x1) - (yt - y1)
    # P2 == -P1 or visa versa
    return xt - x1


"""
Miller's algorithm for elliptic curves
 - Section 2.3, Algorithm 1 (pg 5)
 - https://eprint.iacr.org/2008/096.pdf

Pairings for beginners, Craig Costello
 - Section 5.3, Algorithm 5.1 (pg 79)
 - http://www.craigcostello.com.au/pairings/PairingsForBeginners.pdf

Compute the function f_{s,P} for s > 0
"""


def ate_miller_loop(Q, P, ate_loop_count, field_class):
    """
    Ate pairing
    a_T : G_2 x G_1 -> G_3
    (Q,P) -> f_{T,Q}(P)^{p^k-1/r}
    """
    if not Q or not P:
        return field_class.one()
    assert P != Q
    R = Q
    f = field_class.one()
    # Loop is executed for all bits (except the MSB itself),
    # in MSB to LSB order, skipping all leading zeros
    log_ate_loop_count = floor(log2(ate_loop_count)) - 1
    for i in range(log_ate_loop_count, -1, -1):
        f = (f*f) * linefunc(R, R, P)       # doubling step
        R = R.double()
        if ate_loop_count & (2**i):         # addition step
            f = f * linefunc(R, Q, P)
            R = R.add(Q)
    return f


def ate_pairing(Q: AbstractPointG2, P: AbstractPointG1, group: AbstractGroup, ate_loop_count: int):
    field_class = group.GT().field()
    field_modulus = field_class.field_modulus
    curve_order = group.order()
    degree = field_class.degree
    assert degree % 2 == 0  # degree must be even
    pkm1 = (field_modulus**degree) - 1
    assert (pkm1 % curve_order) == 0

    # Final exponentiation after miller loop
    return ate_miller_loop(Q, P, ate_loop_count, field_class)**(pkm1//curve_order)

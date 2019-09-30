from math import log2, floor

from ..group import AbstractGroup, AbstractPoint, AbstractPointG2, AbstractPointG1
from ..field import make_Fq, make_Fqk
from ..sw import ShortWeierstrassPoint

"""
A Family of Implementation-Friendly BN Elliptic Curves
 - https://eprint.iacr.org/2010/429.pdf
"""


modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583


Fq = make_Fq(modulus)
Fq2 = make_Fqk(modulus, [1, 0])
Fq12 = make_Fqk(modulus, [82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0])


class ALTBN_254_G1(AbstractPointG1, ShortWeierstrassPoint):
    PARAM_A = Fq.zero()
    PARAM_B = Fq(3)

    @classmethod
    def field(cls):
        return Fq

    @classmethod
    def group(cls):
        return ALTBN_254

    @classmethod
    def generator(cls):
        return cls(1, 2)


class ALTBN_254_G2(AbstractPointG2, ShortWeierstrassPoint):
    PARAM_A = Fq2.zero()
    PARAM_B = Fq2([3, 0]) / Fq2([9, 1])

    @classmethod
    def field(cls):
        return Fq2

    @classmethod
    def group(cls):
        return ALTBN_254

    @classmethod
    def generator(cls):
        x = Fq2([
            10857046999023057135944570762232829481370756359578518086990519993285655852781,
            11559732032986387107991004021392285783925812861821192530917403151452391805634,
        ])
        y = Fq2([
            8495653923123431417604973247489272438418190587263600148770280649306958101930,
            4082367875863433681332203403145435568316851327593401208105741076214120093531,
        ])
        return cls(x, y)


class ALTBN_254_GT(ShortWeierstrassPoint):
    PARAM_A = Fq12.zero()
    PARAM_B = Fq12([3] + [0] * 11)

    @classmethod
    def field(cls):
        return Fq12

    @classmethod
    def group(cls):
        return ALTBN_254


class ALTBN_254(AbstractGroup):
    @classmethod
    def order(cls):
        return 21888242871839275222246405745257275088548364400416034343698204186575808495617

    @classmethod
    def G1(cls):
        return ALTBN_254_G1

    @classmethod
    def G2(cls):
        return ALTBN_254_G2

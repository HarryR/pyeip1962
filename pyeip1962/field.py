from typing import Sequence, List, Generator

from math import gcd, log2, ceil

from py_ecc import fields
from py_ecc.fields.field_elements import IntOrFQ

from .redc import mont_findR, mont_convert, mont_redux


def frobenius_coeffs_powers(modulus: int, degree: int, num: int=None, divisor: int=None) -> Generator[int,None,None]:
    divisor = divisor or degree
    num = num or 1
    tower_modulus = modulus ** degree
    for i in range(num):
        a = i + 1
        q_power = 1
        powers = []
        for j in range(degree):
            powers.append((((a*q_power)-a) // divisor) % tower_modulus)
            q_power *= modulus
        yield powers


def frobenius_coeffs(non_residue: int, *args, **kwa) -> List[List[int]]:
    coeffs = list()
    for i, powers in enumerate(frobenius_coeffs_powers(*args, **kwa)):
        coeffs.append(list())
        for p_i in powers:
            coeffs[i].append( non_residue ** p_i )
    return coeffs


def from_hexlimbs(limbs, limb_bits: int = 64) -> int:
    return from_limbs([int(_, 16) for _ in limbs], limb_bits)


def to_limbs(x: int, x_bits: int, limb_bits: int = 64) -> List[int]:
    x = int(x)
    n_limbs = int(ceil(x_bits/limb_bits))
    limb_mask = (1<<limb_bits)-1
    result = list()
    for i in range(n_limbs):
        limb = (x >> (i*limb_bits)) & limb_mask
        result.append(limb)
    return result


def to_hexlimbs(x: int, x_bits: int, limb_bits: int = 64) -> List[str]:
    return [hex(_) for _ in to_limbs(x, x_bits, limb_bits)]


def from_limbs(limbs, limb_bits: int = 64) -> int:
    x = 0
    for i, l_i in enumerate(limbs[::-1]):
        x = x << limb_bits
        x += l_i
    return x


class CommonFieldStuff:
    def norm(self):
        raise NotImplementedError

    def legendre(self):
        euler = (self.field_modulus-1)//2
        s = self.norm() ** euler
        if s == 0:
            return 0
        elif s == 1:
            return 1
        return -1

    def is_quadratic_residue(self):
        return self.legendre() != -1


class AbstractField(fields.FQ, CommonFieldStuff):
    def norm(self):
        return self


class AbstractExtensionField(fields.FQP, CommonFieldStuff):
    def norm(self):
        """
        From: https://eprint.iacr.org/2010/429.pdf '2 Preliminaries':
        The conjugates of `a \in F_{p^e}` are the elements `a^{p^i}`, where `0 <= i < e`.
        The norm N(a) of `a \in F_{p^e}` is defined to be the product of all its conjugates.

            `N(a) := \prod_i a^{p^i}`

        The field norm maps the element onto a smaller ring of its tower.

        References:
         - https://en.wikipedia.org/wiki/Field_norm
         - https://eprint.iacr.org/2019/015.pdf
         - https://file.scirp.org/pdf/JIS_2019070114405856.pdf

        For FP2:
            (c0^2 - beta) * c1^2

        The trace function maps an element of the extension field `F_{p^m}` to
        an element of the prime field `F_p.`
        """
        pass


def make_Fq(q: int):
    R = mont_findR(q)
    q_bits = ceil(log2(q))

    class Fq(AbstractField):
        """Prime field"""
        field_modulus = q
        def __init__(self, val: IntOrFQ) -> None:
            super().__init__(val)

        def __getitem__(self, idx):
            if idx == 0:
                return self.n
            raise KeyError()        

        @classmethod
        def from_limbs(cls, limbs: Sequence[int], limb_bits: int = 64) -> 'Fq':
            return cls(from_limbs(limbs, limb_bits))

        @classmethod
        def from_mont_limbs(cls, limbs: Sequence[int], limb_bits: int = 64) -> 'Fq':
            x = from_limbs(limbs, limb_bits)
            return cls(mont_redux(x, q, R))

        def to_limbs(self, limb_bits: int = 64) -> List[int]:
            return to_limbs(self.n, limb_bits)

        def to_mont_limbs(self, limb_bits: int = 64) -> List[int]:
            return to_limbs(mont_convert(self.n, q, R), q_bits)

    return Fq


def make_Fqk(field_modulus: int, modulus_coeffs: Sequence[IntOrFQ]):
    R = mont_findR(field_modulus)
    q_bits = ceil(log2(field_modulus))
    q = field_modulus

    class Fqk(AbstractExtensionField):
        """Extension field of F_q^k"""
        degree = len(modulus_coeffs)
        field_modulus = q

        def __init__(self,
                     coeffs: Sequence[IntOrFQ]) -> None:
            if isinstance(coeffs, Fqk):
                # Shallow copy
                assert len(coeffs.coeffs) == len(modulus_coeffs)
                assert coeffs.field_modulus == field_modulus
                coeffs = coeffs.coeffs
            super().__init__(coeffs, modulus_coeffs)

        def __getitem__(self, idx):
            if idx < len(self.coeffs):
                return self.coeffs[idx]
            raise KeyError()

        @classmethod
        def from_limbs(cls, limbs: Sequence[Sequence[int]], limb_bits: int = 64) -> 'Fqk':
            return cls([from_limbs(_, limb_bits) for _ in limbs])

        @classmethod
        def from_mont_limbs(cls, limbs: Sequence[Sequence[int]], limb_bits: int = 64) -> 'Fqk':
            mont_coeffs = [from_limbs(_, limb_bits) for _ in limbs]
            return cls([mont_redux(_, field_modulus, R) for _ in mont_coeffs])

        def to_limbs(self, limb_bits: int = 64) -> List[int]:
            return [to_limbs(_, limb_bits) for _ in self.coeffs]

        def to_mont_limbs(self, limb_bits: int = 64) -> List[int]:
            x = [mont_convert(_, field_modulus, R) for _ in self.coeffs]
            return [to_limbs(_, q_bits) for _ in x]

    return Fqk


def test_bls12_377():
    q = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
    q_bits = 377
    R = mont_findR(q)
    Fq = make_Fq(q)
    fp_non_residue = Fq(-5)
    fp2_coeffs = frobenius_coeffs(fp_non_residue, q, 2)[0]
    fp2_coeffs_mont = [to_limbs(mont_convert(_, q, R), q_bits) for _ in fp2_coeffs]
    expected_c1 = [
        [0x2cdffffffffff68,
        0x51409f837fffffb1,
        0x9f7db3a98a7d3ff2,
        0x7b4e97b76e7c6305,
        0x4cf495bf803c84e8,
        0x8d6661e2fdf49a],
        [0x823ac00000000099,
        0xc5cabdc0b000004f,
        0x7f75ae862f8c080d,
        0x9ed4423b9278b089,
        0x79467000ec64c452,
        0x120d3e434c71c50]]
    assert expected_c1 == fp2_coeffs_mont

    # Quadratic non-residue
    qnr_mont = [0] + [from_limbs([
            202099033278250856,
            5854854902718660529,
            11492539364873682930,
            8885205928937022213,
            5545221690922665192,
            39800542322357402])]
    qnr = [mont_redux(_, q, R)[0] for _ in qnr_mont]


def test_bls12_381():
    q = from_limbs([
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a])
    q_bits = ceil(log2(q))
    assert q_bits == 381

    R = mont_findR(q)
    Fq = make_Fq(q)
    non_residue = Fq(-1)
    fp2_coeffs = frobenius_coeffs(non_residue, q, 2)[0]
    fp2_coeffs_mont = [_.to_mont_limbs() for _ in fp2_coeffs]
    expected_c1 = [
        [0x760900000002fffd,
         0xebf4000bc40c0002,
         0x5f48985753c758ba,
         0x77ce585370525745,
         0x5c071a97a256ec6d,
         0x15f65ec3fa80e493],
        [0x43f5fffffffcaaae,
         0x32b7fff2ed47fffd,
         0x7e83a49a2e99d69,
         0xeca8f3318332bb7a,
         0xef148d1ea0f4c069,
         0x40ab3263eff0206]]
    assert fp2_coeffs_mont == expected_c1

    fp2_modulus_coeffs = [1,0]
    Fq2 = make_Fqk(q, fp2_modulus_coeffs)

    # Quadratic non-residue
    qnr_mont = [
        from_limbs([
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493]),
        from_limbs([
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493])]
    qnr = [mont_redux(_, q, R)[0] for _ in qnr_mont]
    assert qnr == [1,1]
    qnr = Fq2(qnr)


    fq6_coeffs_raw = frobenius_coeffs(qnr, q, 6, 2, 3)
    print('fq6_coeffs_raw', fq6_coeffs_raw)
    fq6_coeffs = [[x.to_mont_limbs() for x in _]
                  for _ in fq6_coeffs_raw]
    print('fq6_coeffs', fq6_coeffs)
    for a in fq6_coeffs:
        for b in a:
            for c in b:
                print([hex(_) for _ in c])
    assert fq6_coeffs[0][0] == [
        [0x760900000002fffd,
         0xebf4000bc40c0002,
         0x5f48985753c758ba,
         0x77ce585370525745,
         0x5c071a97a256ec6d,
         0x15f65ec3fa80e493],
        [0, 0, 0, 0, 0, 0]]
    assert fq6_coeffs[1][2] == [
        [0xcd03c9e48671f071,
         0x5dab22461fcda5d2,
         0x587042afd3851b95,
         0x8eb60ebe01bacb9e,
         0x3f97d6e83d050d2,
         0x18f0206554638741],
        [0, 0, 0, 0, 0, 0]]

    fq12_coeffs_raw = frobenius_coeffs(qnr, q, 12, 1, 6)
    print('fq6_coeffs_raw', fq12_coeffs_raw)
    fq12_coeffs = [[x.to_mont_limbs() for x in _]
                  for _ in fq12_coeffs_raw]
    print('fq6_coeffs', fq12_coeffs)
    for a in fq12_coeffs:
        for b in a:
            for c in b:
                print([hex(_) for _ in c])


if __name__ == "__main__":
    from math import ceil

    test_bls12_377()
    test_bls12_381()

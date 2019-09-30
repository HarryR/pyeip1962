from typing import NamedTuple, Tuple, List
import enum
from .fields import Fq, Fqk


@enum.unique
class Operation(enum.IntEnum):
	G1_ADD = 1
	G1_MUL = 2
	G1_MULTIEXP = 3
	G2_ADD = 4
	G2_MUL = 5
	G2_MULTIEXP = 6
	PAIRING = 7

	@property
	def is_g1(self):
		return self.value in [self.G1_ADD, self.G1_MUL, self.G1_MULTIEXP]

	@property
	def is_g2(self):
		return self.value in [self.G2_ADD, self.G2_MUL, self.G2_MULTIEXP]

	@property
	def is_pairing(self):
		return self.value == self.PAIRING


@enum.unique
class CurveFamily(enum.IntEnum):
	BLS12 = 1
	BN = 2
	MNT4 = 3
	MNT6 = 4
	CP = 5


@enum.unique
class TwistType(enum.IntEnum):
	M = 1
	D = 2


class PointG1(NamedTuple):
	x: Fq
	y: Fq


class PointG2(NamedTuple):
	x: Fqk
	y: Fqk


Point = Union[PointG1,PointG2]



class G1Prefix(NamedTuple):
	field_length: int
	field_modulus: int
	A: int
	B: int
	order_length: int
	order: int


class G1AddOp(NamedTuple):
	prefix: G1Prefix
	lhs: G1Point
	rhs: G1Point


class G1MulOp(NamedTuple):
	prefix: G1Prefix
	lhs: G1Point
	rhs: int

class G1MultiExpOp(NamedTuple):
	prefix: G1Prefix
	num_pairs: int
	points: List[Tuple[int,Point]]
	scalar: int


G1Op = Union[G1AddOp, G1MulOp, G1MultiExpOp]


class G2Prefix(NamedTuple):
	field_length: int
	field_modulus: int
	extension_degree: int
	fp_non_residue: int
	A: int
	B: int
	order_length: int
	order: int


class G2AddOp(NamedTuple):
	prefix: G2Prefix
	lhs: G2Point
	rhs: G2Point


class G2MulOp(NamedTuple):
	prefix: G2Prefix
	lhs: G2Point
	scalar: int


class G2MultiExp(NamedTuple):
	prefix: G2Prefix
	num_pairs: int
	point: G2Point
	scalar: List[int]


G2Op = Union[G2AddOp, G2MulOp, G2MultiExp]


class PairingOp(NamedTuple):
	curve_type: CurveFamily
	field_length: int
	field_modulus: int
	A: int
	B: int
	order_length: int
	order: int
	fp2_non_residue: int
	fp6_non_residue: int
	twist_type: TwistType
	x_length: int
	x: int
	sign: int
	num_pairs: int
	pairs: List[Tuple[G1Point,G2Point]]


AnyOp = Union[G1Op, G2Op, PairingOp]


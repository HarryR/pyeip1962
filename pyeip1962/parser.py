from typing import Union, Callable, Any

from .ops import AnyOp, Operation, G1Prefix, G2Prefix, G1AddOp, G1MulOp, G1MultiExpOp, G1Op
from .ops import G2AddOp, G2MulOp, G2MultiExpOp, G2Op
from .ops import G1Point, G2Point, TwistType, CurveFamily
from .structs import 
from .fields import Fq, Fqk


def _make_int(data: bytes) -> int:
	return int.from_bytes(data, 'big')


CURVE_TYPE_LENGTH = 1
OPERATION_ENCODING_LENGTH = 1
TWIST_TYPE_LENGTH = 1
EXTENSION_DEGREE_ENCODING_LENGTH = 1
BYTES_FOR_LENGTH_ENCODING = 1


def is_non_nth_root(divisor, modulus):
	if x == 0:
		return False
	quotient, remainder = divmod(modulus-1, divisor)
	if remainder == 0:
		return False
	l = pow(divisor, quotient, modulus)
	return l != 1


class StreamParser(object):
	def __init__(self, data: bytes):
		self.data = data

	def _consume(self, n: int, parsefn: Callable) -> Any:
		if len(calldata) < n:
			raise RuntimeError(f'Could not read {n} bytes, only {len(calldata)} remaining')
		obj = callfn(calldata[:n])
		self.data = calldata[n:]
		return obj

	def _consume_int(self, n_bytes: int) -> int:
		return self._consume(n_bytes, _make_int)

	def _consume_many_int(self, n_bytes: int, count: int) -> List[int]:
		return [self._consume_int(n_bytes) for _ in range(count)]

	def _consume_Fq(self, n_bytes: int, field_modulus: int) -> Fq:
		return Fq(self._consume_int(n_bytes), field_modulus)

	def _consume_Fqk(self, n_bytes: int, extension_degree: int, field_modulus: int, modulus_coeffs) -> Fqk:
		coeffs = [self._consume_int(n_bytes) for _ in range(extension_degree)]
		return Fqk(coeffs, modulus_coeffs, field_modulus)

	def g1_point(self, prefix: G1Prefix) -> G1Point:
		x = self._consume_Fq(prefix.field_length, prefix.field_modulus)
		y = self._consume_Fq(prefix.field_length, prefix.field_modulus)
		return G1Point(x, y)

	def g1_prefix(self) -> G1Prefix:
		field_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		field_modulus, A, B = self._consume_many_int(field_length, 3)
		order_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		order = self._consume_int(order_length)
		return G1Prefix(field_length, field_modulus, A, B, order_length, order)

	def g1_add_op(self, prefix: G1Prefix) -> G1AddOp:
		return G1AddOp(self.g1_point(prefix), self.g1_point(prefix))

	def g1_point_and_scalar(self, prefix: G1Prefix) -> Tuple[G1Point, int]:
		point = self.g1_point(prefix)
		scalar = self._consume_Fq(prefix.order_length, prefix.field_modulus)
		return (point, scalar)

	def g1_mul_op(self, prefix: G1Prefix) -> G1MulOp:
		point, scalar = self.g1_point_and_scalar(prefix)
		return G1MulOp(prefix, point, scalar)

	def g1_multiexp_op(self, prefix: G1Prefix) -> G1MultiExpOp:
		num_pairs = self._consume_int(1)
		pairs = [self.g1_point_and_scalar(prefix) for _ in range(num_pairs)]
		return G1MultiExpOp(prefix, num_pairs, pairs)

	def g1_op(self, op: Operation) -> G1Op:
		prefix = self.g1_prefix()
		if op == Operation.G1_ADD:
			return self.g1_add_op(prefix)
		elif op == Operation.G1_MUL:
			return self.g1_mul_op(prefix)
		elif op == Operation.G1_MULTIEXP:
			return self.g1_multiexp_op(prefix)

	def g2_prefix(self) -> G2Prefix:
		field_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		field_modulus = self._consume_int(field_length)
		extension_degree = self._consume_int(EXTENSION_DEGREE_ENCODING_LENGTH)
		non_residue, A, B = self._consume_many_int(field_length, 3)
		order_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		order = self._consume_int(order_length)
		return G2Prefix(field_length, field_modulus, non_residue, A, B, order_length, order)

	def g2_point(self, prefix: G2Prefix) -> G2Point:
		x_coeffs = self._consume_many_int(prefix.field_length, prefix.extension_degree)
		y_coeffs = self._consume_many_int(prefix.field_length, prefix.extension_degree)
		return G2Point(x_coeffs, y_coeffs)

	def g2_add_op(self, prefix: G2Prefix):
		return G2AddOp(prefix, self.g2_point(prefix), self.g2_point(prefix))

	def g2_point_and_scalar(self, prefix: G2Prefix) -> Tuple[G2Point, int]:
		point = self.g2_point(prefix)
		scalar = self._consume_int(prefix.order_length)
		return (point, scalar)

	def g2_mul_op(self, prefix: G2Prefix) -> G2MulOp:
		return G2MulOp(*self.g2_point_and_scalar(prefix))

	def g2_multiexp_op(self, prefix: G2Prefix) -> G2MultiExpOp:
		num_pairs = self._consume_int(1)
		pairs = [self.g2_point_and_scalar(prefix) for _ in range(num_pairs)]
		return G2MultiExpOp(prefix, num_pairs, pairs)

	def g2_op(self, op: Operation) -> G2Op:
		prefix = self.g2_prefix()
		if op == Operation.G2_ADD:
			return self.g2_add_op(prefix)
		elif op == Operation.G2_MUL:
			return self.g2_mul_op(prefix)
		elif op == Operation.G2_MULTIEXP:
			return self.g2_multiexp_op(prefix)

	def g1_g2_pair(self, g1_prefix: G1Prefix, g2_prefix: G2Prefix):
		pass

	def pairing_bls12(self, field_length: int, field_modulus: int):
		A, B = self._consume_many_int(field_length, 2)
		if A != 0:
			raise ValueError("A parameter must be zero for BLS12 curve")
		order_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		order = self._consume_int(order_length)
		if order == 0:
			raise ValueError("Group order is zero")

		fp_non_residue = self._consume_int(field_length)
		is_not_a_square = is_non_nth_root(fp_non_residue, field_modulus)
		if not is_not_a_square:
			raise ValueError("Non-residue for Fp2 is actually a residue")


	def pairing_op(self, op: Operation) -> PairingOp:
		curve_type = CurveFamily(self._consume_int(CURVE_TYPE_LENGTH))
		field_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		field_modulus = self._consume_int(field_length)

		if curve_type == CurveFamily.BLS12:
			return self.pairing_bls12(field_length, field_modulus)
		elif curve_type == CurveFamily.BN:
			return self.pairing_bn(field_length, field_modulus)
		elif curve_type == CurveFamily.MNT4:
			return self.pairing_mnt4(field_length, field_modulus)
		elif curve_type == CurveFamily.MNT6:
			return self.pairing_mnt6(field_length, field_modulus)

		A, B = self._consume_many_int(field_length, 2)
		order_length = self._consume_int(BYTES_FOR_LENGTH_ENCODING)
		order = self._consume_int(order_length)
		fp2_non_residue = self._consume_int(field_length)
		fp6_non_residue = self._consume_many_int(field_length, 2)
		twist_type = TwistType(self._consume_int(TWIST_TYPE_LENGTH))
		x_length = self._consume_int(1)
		x = self._consume_int(x_length)
		sign, num_pairs = self._consume_many_int(1, 2)
		pairs = 

	def parse(self) -> AnyOp:
		op_code = self._consume(OPERATION_ENCODING_LENGTH, lambda data: Operation(_make_int(data)))
		if op_code.is_g1:
			return self.g1_op(op_code)
		elif op_code.is_g1:
			return self.g2_op(op_code)
		elif op_code.is_pairing:
			return self.pairing_op(op_code)


import unittest


def group_law_tests(group):
	g = group.generator()
	assert g.is_on_curve()

	# Negation and zero works
	assert g is not None
	assert g - g == group.zero()
	assert not (g - g)

	# Order of the curve
	assert g * group.order() == group.zero()
	assert g * (group.order() - 1) == -g

	# Verify doubling works as expected
	g4 = g + g + g + g
	assert g4 == g.double().double()


def pairing_tests(curve):
	g1 = curve.G1().generator()
	g2 = curve.G2().generator()
	g1_20 = g1 * 20
	g2_20 = g2 * 20

	a = curve.pairing(g1_20, g2)
	b = curve.pairing(g1, g2_20)
	assert a == b

	# TODO: more tests for bilinearaity


class GroupTests(unittest.TestCase):
	"""
	def test_altbn_254(self):
		from pyeip1962.curves.altbn_254 import ALTBN_254
		group_law_tests(ALTBN_254.G1())
		group_law_tests(ALTBN_254.G2())
		#pairing_tests(ALTBN_254)

	def test_bls12_381(self):
		from pyeip1962.curves.bls12_381 import BLS12_381
		group_law_tests(BLS12_381.G1())
		group_law_tests(BLS12_381.G2())
		pairing_tests(BLS12_381)
	"""

	def test_bls12_377(self):
		from pyeip1962.curves.bls12_377 import BLS12_377
		group_law_tests(BLS12_377.G1())
		group_law_tests(BLS12_377.G2())
		pairing_tests(BLS12_377)


if __name__ == "__main__":
	unittest.main()

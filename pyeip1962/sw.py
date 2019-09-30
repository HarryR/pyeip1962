from .group import AbstractPoint


class ShortWeierstrassPoint(AbstractPoint):
	"""
	y^2 = x^3 + a4*x + a6
	"""

	PARAM_A = None		# a4
	PARAM_B = None		# a6

	def neg(self):
		return type(self)(self.x, -self.y)

	def add(self, other):
		if not other:
			return self
		elif all([self.x == other.x, self.y == other.y]):
			return self.double()
		elif self.x == other.x:
			# Add self, to its negative, equals zero
			return self.zero()
		# TODO: avoid inversions
		lam = (other.y - self.y) / (other.x - self.x)
		x3 = lam**2 - self.x - other.x
		y3 = (lam*(self.x - x3)) - self.y
		return type(self)(x3, y3)

	def double(self):
		# TODO: avoid inversions
		lam = (3*(self.x**2)+self.PARAM_A) / (2*self.y)
		x3 = lam**2 - (2*self.x)
		y3 = lam*(self.x - x3) - self.y
		return type(self)(x3, y3)

	def is_on_curve(self):
		# y^2=x^3+a*x+b
		ysq = (self.x**3) + (self.PARAM_A*self.x) + self.PARAM_B
		return (self.y**2) == ysq

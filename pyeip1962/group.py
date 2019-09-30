
class AbstractPoint(object):
	__slots__ = ('x', 'y')

	def __init__(self, x, y):
		self.x = self.field()(x)
		self.y = self.field()(y)

	def __neg__(self):
		return self.neg()

	def __add__(self, other):
		return self.add(other)

	def __sub__(self, other):
		return self.add(other.neg())

	def __mul__(self, n):
		return self.mul(n)

	def __iter__(self):
		return iter([self.x, self.y])

	def __eq__(self, other: 'AbstractPoint'):
		return other is not None and self.x == other.x and self.y == other.y

	def __bool__(self):
		return self != self.zero()

	def __repr__(self):
		return f'{type(self).__name__}(x={self.x}, y={self.y})'

	@classmethod
	def generator(cls):
		raise NotImplementedError

	@classmethod
	def zero(cls):
		return None

	@classmethod
	def field(cls):
		# Field used for X and Y coordinates
		raise NotImplementedError

	@classmethod
	def order(cls):
		return cls.group().order()

	@classmethod
	def group(cls):
		# Group
		raise NotImplementedError

	def neg(self):
		raise NotImplementedError

	def add(self, other: 'AbstractPoint'):
		raise NotImplementedError

	def double(self, other: 'AbstractPoint'):
		raise NotImplementedError

	def mul(self, scalar):
		scalar = int(scalar)
		if scalar == 1:
			return self
		p = self
		a = self.zero()
		while scalar != 0:
			if (scalar & 1) != 0:
				a = p.add(a)
			p = p.double()
			scalar = scalar // 2
		return a


class AbstractPointG1(AbstractPoint):
	def pairing(self, other: 'AbstractPointG2'):
		return self.group().pairing(self, other)


class AbstractPointG2(AbstractPoint):
	def pairing(self, other: AbstractPointG1):
		return self.group().pairing(other, self)


class AbstractGroup(object):
	@classmethod
	def order(cls):
		raise NotImplementedError

	@classmethod
	def G1(cls):
		"""Returns class for G1 group"""
		raise NotImplementedError

	@classmethod
	def G2(cls):
		"""Returns class for G2 group"""
		raise NotImplementedError

	@classmethod
	def GT(cls):
		"""Returns class for target group"""
		raise NotImplementedError

	@classmethod
	def pairing(cls, a: AbstractPointG1, b: AbstractPointG2):
		assert isinstance(a, self.G1())
		assert isinstance(b, self.G2())
		raise NotImplementedError

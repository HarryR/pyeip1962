# Computes the GCD of a and b.
def gcd(a, b):
	while b:
		a, b = b, a%b
	return a

# Computes XGCD(a, b).
# Returns (m, x, y) s.t. m = a*x + b*y
def xgcd(a, b):
	prev_x, x = 1, 0
	prev_y, y = 0, 1
	while b:
		q, a, b = a//b, b, a % b
		x, prev_x = prev_x - q*x, x
		y, prev_y = prev_y - q*y, y
	return a, prev_x, prev_y

# Finds an R such that R = 2^k, R > N, for the smallest k.
def mont_findR(N, limb_size=64):
	g = 0
	b = 2 ** limb_size
	R = b
	while g != 1:
		R *= b
		if R > N:
			g = gcd(R, N)
	return R

# Converts T into Montgomery Form.
def mont_convert(T, N, R):
	return int(T*R) % N

# Produces the Montgomery Reduction of T modulo N.
# Returns tuple of the reduction and whether or not we reduced mod N
def mont_redux(T, N, R, Ni = None):
	if Ni == None:
		_, _, Ni = xgcd(R, N)
	m = -(T*Ni) % R
	t = ((T + (m*N)) // R)
	return t % N, t > N

# Adds two numbers, a and b, in Montgomery form, under modulo N.
# Returns an answer in Montgomery Form.
def mont_add(a, b, N):
	return int(a + b) % N

# Subtracts two numbers, b from a, in Montgomery form, under modulo N.
# Returns an answer in Montgomery Form.
def mont_sub(a, b, N):
	return int(a - b) % N

# Multiplies two numbers, a and b, in Montgomery form, under modulo N.
# Returns an answer in Montgomery Form in a tuple, with whether or not we
# reduced the answer modulo N.
def mont_mul(a, b, N, R, Ni = None):
	return mont_redux(a * b, N, R, Ni)

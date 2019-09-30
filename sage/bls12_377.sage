field_modulus = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
desired_curve_order = 8444461749428370424248824938781546531375899335154063827935233455917409239041
curve_name = 'bls12_377'

Fp = GF(field_modulus)

PARAM_A4 = Fp(0)
PARAM_A6 = Fp(1)

E = EllipticCurve(Fp, [PARAM_A4, PARAM_A6])
E_order = E.order()
print '# E order', E_order
assert E_order % desired_curve_order == 0
assert desired_curve_order.is_prime() is True
E_cofactor = E_order // desired_curve_order
Fr = GF(desired_curve_order)

# 12th degree extension field must be pairing friendly
assert ((field_modulus**12)-1)%desired_curve_order == 0 

R.<T> = PolynomialRing(Fp)

# Starting at -1 and subtracting is an arbitrary choice
# could start at 1, where 2 will usually be the first non-residue
for non_residue in range(field_modulus-1, field_modulus-(10**3), -1):
    non_residue = Fp(non_residue)
    if non_residue.is_square():
        continue
    F2_equation = T^2-non_residue
    F2.<u> = Fp.extension(F2_equation, 'u')
    success = False
    for j in range(0, 10**3):
        quadratic_non_residue = u+j
        # Must not be square residue
        if quadratic_non_residue.is_square():
            continue
        try:
            # Must not be cubic residue
            quadratic_non_residue.nth_root(3)
            continue
        except ValueError as ex:
            pass
        F12_equation = (T^6 - j)^2 - non_residue
        u_to_w = T^6 - j
        w_to_u = T + j
        success = True
        break
    if success:
        break


print '# non residue', non_residue
print '# quadratic_non_residue', quadratic_non_residue

F12.<w> = Fp.extension(F12_equation)

print '# w', F12.vector_space()(w)
print '# u to w', F12.vector_space()(u_to_w(w))

E12 = EllipticCurve(F12, [PARAM_A4, PARAM_A6])

E2 = EllipticCurve(F2, [PARAM_A4*quadratic_non_residue, PARAM_A6*quadratic_non_residue])
is_D_type = False
A_twist = 0
if not (E2.order()/desired_curve_order).is_integer():
    # E2 requires a twist to be divisible by the desired curve order
    B_twist = PARAM_A6/quadratic_non_residue
    E2 = EllipticCurve(F2, [0, B_twist])
    if not (E2.order()/desired_curve_order).is_integer():
        raise Exception("twist didn't have appropriate order")
    is_D_type = True
    print("# D type twist")
    F2_PARAM_A4 = PARAM_A4 / quadratic_non_residue
    F2_PARAM_A6 = PARAM_A6 / quadratic_non_residue
else:
    # E2 order is exactly divisible by curve order
    B_twist = PARAM_A6*quadratic_non_residue
    F2_PARAM_A6 = PARAM_A6 * quadratic_non_residue
    F2_PARAM_A4 = PARAM_A4 * quadratic_non_residue
    print('# M type twist')


E2_order = E2.order()
assert E2_order % desired_curve_order == 0
E2_cofactor = E2_order // desired_curve_order


def frobenius_coeffs_powers(modulus, degree, num=None, divisor=None):
    divisor = divisor or degree
    num = num or 1
    tower_modulus = modulus ** degree
    for i in range(num):
        a = i + 1
        q_power = 1
        powers = []
        for j in range(degree):
            x = (((a*q_power)-a) // divisor) % tower_modulus
            powers.append(x)
            q_power *= modulus
        yield powers


def frobenius_coeffs(non_residue, modulus, *args, **kwa):
    coeffs = list()
    for i, powers in enumerate(frobenius_coeffs_powers(modulus, *args, **kwa)):
        coeffs.append(list())
        for p_i in powers:
            coeffs[i].append( pow(non_residue, p_i, modulus) )
    return coeffs


def find_generator(E, F, a6, cofactor, order):
    for x in range(1, 10**3):
        x = F(x)
        y2 = x**3 + a6
        if not y2.is_square():
            continue
        y = y2.sqrt()
        p = cofactor*E(x, y)
        if not p.is_zero() and (order*p).is_zero():
            # Choose the 'positive' point, where the Y coordinate is the lowest
            # This is an arbitrary choice... but allows for consistency across implementations
            negy = -p[1]
            if negy < p[1]:
                return -p
            return p


def find_s_t(name, n):
    for s in range(1, 50):
        if n % (2**s) == 0:
            t = n / 2**s
            assert t.is_integer()
            if not ((t-1)/2).is_integer():
                continue
            print name, "s", s
            print name, "t", t
            print name, "t_minus_1_over_2", (t-1)/2
            return s, t


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


def print_field(name, q, F):
    print name, 'modulus', q
    print name, 'num_bits', ceil(log(q,2))
    print name, 'euler', (q-1)//2
    s, t = find_s_t(name, q-1)
    #print '# finding multiplicative_generator'
    #print name, 'multiplicative_generator', F.vector_space()(F.multiplicative_generator())
    gen = F.gen()
    for i in range(0, 100): #for i in range(-1, -100, -1):
        i = gen+i
        if not i.is_square():
            i_to_t = i**t
            print name, 'nqr', F.vector_space()(i)
            print name, 'nqr_to_t', F.vector_space()(i_to_t)
            break
    print ""


def print_R(name, q, nbits):
    R = mont_findR(q, nbits)
    print name, "R (%d-bit)" % (nbits,), R
    print name, "Rcubed (%d-bit)" % (nbits,), (R**2) % q
    print name, "Rsquared (%d-bit)" % (nbits,), (R**3) % q    
    print name, "inv (%d-bit)" % (nbits,), hex((-q^-1) % 2**nbits)


print curve_name, 'G1 A', Fp.vector_space()(PARAM_A4)
print curve_name, 'G1 B', Fp.vector_space()(PARAM_A6)
print curve_name, 'G1 cofactor', E_cofactor
print curve_name, 'G1 cofactor^1 mod r', (E_cofactor^-1) % desired_curve_order
G1 = find_generator(E, Fp, PARAM_A6, E_cofactor, desired_curve_order)
print curve_name, 'G1_zero', E(0)
print curve_name, 'G1_one', [Fp.vector_space()(_) for _ in G1]
print ""

print curve_name, 'G2 A', F2.vector_space()(F2_PARAM_A4)
print curve_name, 'G2 B', F2.vector_space()(F2_PARAM_A6)
print curve_name, 'G2 cofactor', E2_cofactor
print curve_name, 'G2 cofactor^1 mod r', (E2_cofactor^-1) % desired_curve_order
G2 = find_generator(E2, F2, F2_PARAM_A6, E2_cofactor, desired_curve_order)
print curve_name, 'G2_zero', E2(0)
print curve_name, 'G2_one', [F2.vector_space()(_) for _ in G2]
print ""

print_field(' '.join([curve_name, 'Fq2']), field_modulus**2, F2)

#print '# Calculating Fp^2 coeffs'
print '# F2 polynomial coeffs:', F2_equation.coefficients(sparse=False)[:2]
fp2_coeffs = frobenius_coeffs(non_residue, field_modulus, 2)
for i, c in enumerate(fp2_coeffs):
    print curve_name, 'Fq2', 'Frobenius_coeffs_c1[%d]' % (i,), '=', c
print ''

#print '# Calculating Fp^6 coeffs'
fp6_coeffs = frobenius_coeffs(quadratic_non_residue, field_modulus, 6, 2, 3)
for i, _ in enumerate(fp6_coeffs):
    for j, c in enumerate(_):
        print curve_name, 'Fq6', 'Frobenius_coeffs_c%d[%d]' % (i + 1, j), '=', F2.vector_space()(c)
print ''

# Fq12 is two 6th degree towers
#print '# Calculating Fp^12 frobenius coeffs'
print '# F12 polynomial coeffs:', F12_equation.coefficients(sparse=False)[:12]
fp12_coeffs = frobenius_coeffs(quadratic_non_residue, field_modulus, 12, 1, 6)
for i, _ in enumerate(fp12_coeffs):
    for j, c in enumerate(_):
        print curve_name, 'Fq12', 'Frobenius_coeffs_c%d[%d]' % (i + 1, j), '=', F2.vector_space()(c)
print ''

print_R(' '.join([curve_name, 'Fr']), desired_curve_order, 64)
print_R(' '.join([curve_name, 'Fr']), desired_curve_order, 32)
print_field(' '.join([curve_name, 'Fr']), desired_curve_order, Fr)

print_R(' '.join([curve_name, 'Fq']), field_modulus, 64)
print_R(' '.join([curve_name, 'Fq']), field_modulus, 32)
print_field(' '.join([curve_name, 'Fq']), field_modulus, Fp)
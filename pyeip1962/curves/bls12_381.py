from math import log2, floor

from ..group import AbstractGroup, AbstractPointG2, AbstractPointG1
from ..field import make_Fq, make_Fqk
from ..sw import ShortWeierstrassPoint
from ..pairing import ate_pairing


modulus = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787

Fq = make_Fq(modulus)
Fq2 = make_Fqk(modulus, [1, 0])
Fq12 = make_Fqk(modulus, [2, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0])


class BLS12_381_G1(AbstractPointG1, ShortWeierstrassPoint):
    PARAM_A = Fq.zero()
    PARAM_B = Fq(4)

    @classmethod
    def field(cls):
        return Fq

    @classmethod
    def group(self):
        return BLS12_381

    @classmethod
    def generator(cls):
        x = Fq(3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507)
        y = Fq(1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569)
        return cls(x, y)

    def cast_point_to_fq12(self):
        return BLS12_381_GT(Fq12([self.x] + [0] * 11), Fq12([self.y] + [0] * 11))


class BLS12_381_G2(AbstractPointG2, ShortWeierstrassPoint):
    PARAM_A = Fq2.zero()
    PARAM_B = Fq2([4, 4])

    @classmethod
    def field(cls):
        return Fq2

    @classmethod
    def group(self):
        return BLS12_381

    @classmethod
    def generator(cls):
        x = Fq2([
            352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160,
            3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758,
        ])
        y = Fq2([
            1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905,
            927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582,
        ])
        return cls(x, y)

    def twist_to_GT(self):
        # "Twist" a point in E(FQ2) into a point in E(FQ12)
        w = Fq12([0, 1] + [0] * 10)

        # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 2*x + 2
        x0, x1 = [self.x[0] - self.x[1], self.x[1]]
        y0, y1 = [self.y[0] - self.y[1], self.y[1]]
        # Isomorphism into subfield of Z[p] / w**12 - 2 * w**6 + 2,
        # where w**6 = x
        nx = Fq12([x0] + [0] * 5 + [x1] + [0] * 5)
        ny = Fq12([y0] + [0] * 5 + [y1] + [0] * 5)
        # Divide x coord by w**2 and y coord by w**3
        return BLS12_381_GT(nx / w**2, ny / w**3)


class BLS12_381_GT(ShortWeierstrassPoint):
    PARAM_A = Fq12.zero()

    @classmethod
    def field(cls):
        return Fq12

    @classmethod
    def generator(cls):
        return BLS12_381_G2.generator().twist_to_GT()

    @classmethod
    def group(cls):
        return BLS12_381


class BLS12_381(AbstractGroup):
    @classmethod
    def order(cls):
        return 52435875175126190479447740508185965837690552500527637822603658699938581184513

    @classmethod
    def G1(cls):
        return BLS12_381_G1

    @classmethod
    def G2(cls):
        return BLS12_381_G2

    @classmethod
    def GT(cls):
        return BLS12_381_GT

    @classmethod
    def pairing(cls, a: BLS12_381_G1, b: BLS12_381_G2):
        Q = b.twist_to_GT()
        P = a.cast_point_to_fq12()
        ate_loop_count = 15132376222941642752
        return ate_pairing(Q, P, cls, ate_loop_count)

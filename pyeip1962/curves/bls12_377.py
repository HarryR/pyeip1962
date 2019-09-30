from math import log2, floor

from ..group import AbstractGroup, AbstractPointG2, AbstractPointG1
from ..field import make_Fq, make_Fqk
from ..sw import ShortWeierstrassPoint
from ..pairing import ate_pairing


modulus = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177

Fq = make_Fq(modulus)
Fq2 = make_Fqk(modulus, [5, 0])
Fq12 = make_Fqk(modulus, [5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


class BLS12_377_G1(AbstractPointG1, ShortWeierstrassPoint):
    PARAM_A = Fq.zero()
    PARAM_B = Fq(1)

    @classmethod
    def field(cls):
        return Fq

    @classmethod
    def group(self):
        return BLS12_377

    @classmethod
    def generator(cls):
        x = Fq(81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695)
        y = Fq(17397676153253620270863855454307851802466321586312764156125140564607560990561071773762088186709545111705113293147)
        return cls(x, y)

    def cast_point_to_fq12(self):
        return BLS12_377_GT(Fq12([self.x] + [0] * 11), Fq12([self.y] + [0] * 11))


class BLS12_377_G2(AbstractPointG2, ShortWeierstrassPoint):
    PARAM_A = Fq2.zero()
    PARAM_B = Fq2([0, 155198655607781456406391640216936120121836107652948796323930557600032281009004493664981332883744016074664192874906])

    @classmethod
    def field(cls):
        return Fq2

    @classmethod
    def group(self):
        return BLS12_377

    @classmethod
    def generator(cls):
        # Deterministically derived
        x = Fq2([
            39292833563790338514455678255839969442444299076493345799525535236324569704972737101027043002275594504529645125033,
            97668274349181098911216378040700666521757961257997861327997265570326738925466145318868002777904267769221513117576,
        ])
        y = Fq2([
            12670168495311570839246849220246345469108307986667888010668101126790399240749545663887747620979098015764659835358,
            84432745052336413615082002597703423810618940985259643064855840274752478639694687962835382580467718604598437838768,
        ])
        """
        # G2 from Zexe BLS12_377
        x = Fq2([
            233578398248691099356572568220835526895379068987715365179118596935057653620464273615301663571204657964920925606294,
            140913150380207355837477652521042157274541796891053068589147167627541651775299824604154852141315666357241556069118,
        ])
        y = Fq2([
            63160294768292073209381361943935198908131692476676907196754037919244929611450776219210369229519898517858833747423,
            149157405641012693445398062341192467754805999074082136895788947234480009303640899064710353187729182149407503257491,
        ])
        """
        return cls(x, y)

    def twist_to_GT(self):
        # "Twist" a point in E(FQ2) into a point in E(FQ12)
        w = Fq12([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])

        # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 2*x + 2
        x0, x1 = [self.x[0] - self.x[1], self.x[1]]
        y0, y1 = [self.y[0] - self.y[1], self.y[1]]
        # Isomorphism into subfield of Z[p] / w**12 - 2 * w**6 + 2,
        # where w**6 = x
        nx = Fq12([x0] + [0] * 5 + [x1] + [0] * 5)
        ny = Fq12([y0] + [0] * 5 + [y1] + [0] * 5)
        # Divide x coord by w**2 and y coord by w**3
        return BLS12_377_GT(nx * w**2, ny * w**3)


class BLS12_377_GT(ShortWeierstrassPoint):
    PARAM_A = Fq12.zero()

    @classmethod
    def field(cls):
        return Fq12

    @classmethod
    def generator(cls):
        return BLS12_377_G2.generator().twist_to_GT()

    @classmethod
    def group(cls):
        return BLS12_377


class BLS12_377(AbstractGroup):
    @classmethod
    def order(cls):
        return 8444461749428370424248824938781546531375899335154063827935233455917409239041

    @classmethod
    def G1(cls):
        return BLS12_377_G1

    @classmethod
    def G2(cls):
        return BLS12_377_G2

    @classmethod
    def GT(cls):
        return BLS12_377_GT

    @classmethod
    def pairing(cls, a: BLS12_377_G1, b: BLS12_377_G2):
        Q = b.twist_to_GT()
        P = a.cast_point_to_fq12()
        ate_loop_count = 0x8508c00000000001
        return ate_pairing(Q, P, cls, ate_loop_count)

from curve_affine import ECC, Point, O_POINT_INF
import math

# solve the discrete log problem for Q = k*P
def bsgs(curve, P, Q, n=None):

    if not curve.onCurve(Q):
        raise ValueError("Q not on curve: "+str(hex(Q)))
    if not curve.onCurve(P):
        raise ValueError("P not on curve: "+str(hex(P)))
    
    # max, we can supply smaller scope for testing
    if n is None:
        n = curve.n
    sqn = int(math.sqrt(n)) + 1

    # baby steps
    R = O_POINT_INF 
    prec_hash = { O_POINT_INF : 0 }
    for alpha in range(1, sqn):
        R = curve.add(R, P)
        prec_hash[R] = alpha

    # giant steps
    R = Q
    invP = curve.inv(P)
    S = curve.scalarMult(sqn, invP)

    #print(prec_hash)

    for i in range(0, sqn):
        try:
            alpha = prec_hash[R]
        except KeyError:
            pass
        else:
            steps = sqn + i
            lg = alpha + sqn * i
            return lg

        R = curve.add(R, S)

    raise ValueError("failed to solve discrete log problem")

"""
curve = ECC()
P = Point(0xbc90005abc661767049e5653986be754219a67793b4b167a4fbd92e574db8168,
         0x7fa7d9bf70c98557682d6c3686b8ec6b8992aa1e81a47c3257a4c4c527138cf2)
Q = Point(0x60ba5a5aa10ef17557d69e37d11338eea2285a83d8e3a3feb7b15930fe704b8a,
         0x7c2c0cc55a3557d08d045bc14b6657e29a3ef8b7da240934e7b0a58bf5d3c12e)

x = 0xa3b1799e
n = 2**32
print("Curve order: "+str(curve.n))
print("P: "+str(hex(P.x)) + "," + str(hex(P.y)))
print("Q: "+str(hex(Q.x)) + "," + str(hex(Q.y)))
print(str(hex(x)) + "* P = Q")

y = bsgs(curve, P, Q, n)
print('log(p, q) =', y)
print("as hex: "+str(hex(y)))
"""

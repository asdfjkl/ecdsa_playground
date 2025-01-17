from curve_affine import ECC, Point
import random
import sys
from bsgs import bsgs
from tqdm import tqdm

def prepend(s, desired_len):
    prep = ""
    for i in range(0, desired_len - len(s)):
        prep += "0"
    return prep + s

# generate random values k, l (in 32 bit range)
# and compute P = l*G and Q = k * P = k * l * G
# then use baby-step giant-step to search for k

def test_bsgs():
    random.seed(42)
    nist = ECC()
    for i in (range(0,10)): # increase for more test coverage
        k = random.randint(1, 0xFFFFFFFF)
        l = random.randint(0x1, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF)
        P = nist.scalarMult(l, nist.basePoint)
        Q = nist.scalarMult(k, P)

        kl = (k * l) % nist.n

        leak_val_hex = prepend(str(hex(kl))[2:], 64)
        sig_r_hex = prepend(str(hex(P.x))[2:], 64)
        mask_start_hex = prepend(str(hex(k-10))[2:], 8)
        mask_stop_hex = prepend(str(hex(k + 10))[2:], 8)
        #print(leak_val_hex + " "+sig_r_hex + " "+ mask_start_hex + " "+mask_stop_hex)

        #mask_hex = prepend(str(hex(k))[2:], 64)
        #print("mask:        "+str(hex(k)))
        #print("leak_val:    "+ str(hex(kl)))
        #print("r:           "+str(hex(P.x)))
        #print("P.y:         "+str(hex(P.y)))
        #print("Qx:          "+str(hex(Q.x)))
        #print("Qy:          "+str(hex(Q.y)))


        QQ = nist.scalarMult(kl, nist.basePoint)

        assert(Q.x == QQ.x and Q.y == QQ.y)

        # run baby-step giant-step
        # but instead of searchin in the whole range of the 
        # cyclic group (i.e. here nist parameter n) limit
        # to 32 bit range (for test purposes)
        lg = bsgs(nist, P, Q, n=2**32)
        assert(lg == k)


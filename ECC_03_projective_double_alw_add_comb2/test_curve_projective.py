from .curve_projective import Point, ProPoint, ProECC, O_POINT_INF

class TestECCProjective():

    # all tests here are on the default (Nist P-256) curve
    def testAdditions(self):
        ecc = ProECC()

        # bp + O
        proBp = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        pointBp = ecc.aff2pro(Point(0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
                       0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5),42)
        assert(ecc.isEqual(ecc.add(proBp, O_POINT_INF), pointBp))

        # bp + bp
        proBp1 = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        proBp2 = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        point2Bp = ecc.aff2pro(Point(0x7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978,
                       0x07775510DB8ED040293D9AC69F7430DBBA7DADE63CE982299E04B79D227873D1),42)
        assert (ecc.isEqual(ecc.add(proBp1, proBp2), point2Bp))

        # bp + bp + bp
        point1 = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        point2 = point2Bp
        point3Bp = ecc.aff2pro(Point(0x5ECBE4D1A6330A44C8F7EF951D4BF165E6C6B721EFADA985FB41661BC6E7FD6C,
                         0x8734640C4998FF7E374B06CE1A64A2ECD82AB036384FB83D9A79B127A27D5032),42)
        assert (ecc.isEqual(ecc.add(point1, point2), point3Bp))

        # bp + bp + bp + bp
        point1 = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        point2 = point3Bp
        point4Bp = ecc.aff2pro(Point(0xE2534A3532D08FBBA02DDE659EE62BD0031FE2DB785596EF509302446B030852,
                 0xE0F1575A4C633CC719DFEE5FDA862D764EFC96C3F30EE0055C42C23F184ED8C6), 42)
        assert (ecc.isEqual(ecc.add(point1, point2), point4Bp))


    def testScalarMultiplication(self):

        ecc = ProECC()
        # compute 4 * bp by comb2 unmasked
        ecc.proBasePoint = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)

        bp = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        BP4 = ecc.scalarMultComb2Unmasked(4, bp)
        pointRes = ecc.aff2pro(Point(0xE2534A3532D08FBBA02DDE659EE62BD0031FE2DB785596EF509302446B030852,
                         0xE0F1575A4C633CC719DFEE5FDA862D764EFC96C3F30EE0055C42C23F184ED8C6), 42)
        assert(ecc.isEqual(pointRes, BP4))

        ecc = ProECC()
        # compute 4 * bp by comb2 masked
        ecc.proBasePoint = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        bp = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        BP4 = ecc.scalarMultComb2Masked(4, bp, 0)
        pointRes = ecc.aff2pro(Point(0xE2534A3532D08FBBA02DDE659EE62BD0031FE2DB785596EF509302446B030852,
                         0xE0F1575A4C633CC719DFEE5FDA862D764EFC96C3F30EE0055C42C23F184ED8C6), 42)
        assert(ecc.isEqual(pointRes, BP4))


        ecc = ProECC()
        # compute n * bp
        bp = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        BPmul = ecc.scalarMultComb2Masked(112233445566778899, bp, 0)
        pointRes = ecc.aff2pro(Point(0x339150844EC15234807FE862A86BE77977DBFB3AE3D96F4C22795513AEAAB82F,
                         0xB1C14DDFDC8EC1B2583F51E85A5EB3A155840F2034730E9B5ADA38B674336A21), 42)
        assert(ecc.isEqual(pointRes, BPmul))


        ecc = ProECC()
        # compute n * bp
        bp = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
        BPmul = ecc.scalarMultComb2Masked(29852220098221261079183923314599206100666902414330245206392788703677545185283, bp, 0)
        pointRes = ecc.aff2pro(Point(0x9EACE8F4B071E677C5350B02F2BB2B384AAE89D58AA72CA97A170572E0FB222F,
                         0x1BBDAEC2430B09B93F7CB08678636CE12EAAFD58390699B5FD2F6E1188FC2A78), 42)
        assert(ecc.isEqual(pointRes, BPmul))


    def testScalarMultiplicationVecs(self):
        k=None
        x=None
        y=None
        with open("../testvecs_p256.txt", "r") as f:
            for line in f:
                if(k!= None and x!= None and y!= None):
                    # do computation
                    ecc = ProECC()
                    bp = ecc.aff2pro(ecc.pro2aff(ecc.proBasePoint), 42)
                    BPMul = ecc.scalarMultComb2Masked(k, bp, 0)
                    pointRes = ecc.aff2pro(Point(x,y),42)
                    assert (ecc.isEqual(pointRes, BPMul))
                    #print("k: "+str(k))
                    #print("x: "+str(x))
                    #print("y: "+str(y))
                    k = x = y = None
                var_val = line.split("=")
                if (len(var_val) == 2):
                    var = var_val[0].strip()
                    val = var_val[1].strip()
                    if (var == "k"):
                        k = int(val)
                    if(var == "x"):
                        x = int(val, 16)
                    if(var == "y"):
                        y = int(val, 16)

tester = TestECCProjective()
tester.testScalarMultiplicationVecs()
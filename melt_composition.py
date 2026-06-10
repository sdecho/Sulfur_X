# This script decsribes the melt composition change as a function of degree of melt fraction if crystallization is
# enabled (choice == 1). The example below is for Fuego magma. If crystallization is disabled, melt composition does
# not change. The example below is for Hawaiian magma.
# The Please revise this accordingly for the composition of the volcanic system of interest, and notify the changes
# if the results are in use of a publication.

class MeltComposition:
    """melt fraction: 0-1
        choice (==1): crystallization is enabled
    """

    def __init__(self, melt_fraction, choice):

        if choice == 1:# Input melt composition change as a function of k2O here if crystallization is enabled.
             if melt_fraction > 0.16:
                k2o = 0.5 / melt_fraction
                sio2 = -1.9877 * k2o ** 2 +14.864 * k2o + 44.622
                al2o3 = -0.8579 * k2o ** 2 +2.7875 * k2o+ 14.191
                feot = 1.2447 * k2o**2-8.13*k2o + 15.65
                mgo = 1.6993 * k2o ** 2 - 8.9497 * k2o+ 12.026
                cao = 0.1883 * k2o ** 2 - 2.8933 * k2o + 8.8376
                na2o = 5
                p2o5 = -0.0124 * k2o ** 2 - 0.2644 * k2o + 0.7405
                mno = 0.0718 * k2o ** 2-0.3772*k2o + 0.5545
                tio2 = 0.1503 * k2o ** 2 - 0.8544 * k2o + 1.6476
             else:
                melt_fraction = 0.16
                k2o = 0.5 / melt_fraction
                sio2 = -1.9877 * k2o ** 2 + 14.864 * k2o + 44.622
                al2o3 = -0.8579 * k2o ** 2 + 2.7875 * k2o + 14.191
                feot = 1.2447 * k2o ** 2 - 8.13 * k2o + 15.65
                mgo = 1.6993 * k2o ** 2 - 8.9497 * k2o + 12.026
                cao = 0.1883 * k2o ** 2 - 2.8933 * k2o + 8.8376
                na2o = 5
                p2o5 = -0.0124 * k2o ** 2 - 0.2644 * k2o + 0.7405
                mno = 0.0718 * k2o ** 2 - 0.3772 * k2o + 0.5545
                tio2 = 0.1503 * k2o ** 2 - 0.8544 * k2o + 1.6476
            
        else: # Input melt composition here if crystallization is disabled.
            # Fogo_comparison paper
            # sio2 = 42.40 / melt_fraction
            # al2o3 = 11.17 / melt_fraction
            # feot = 12.00 / melt_fraction
            # mgo = 9.55 / melt_fraction
            # cao = 13.31 / melt_fraction
            # na2o = 3.36 / melt_fraction
            # k2o = 1.57 / melt_fraction
            # p2o5 = 0.75 / melt_fraction
            # mno = 0.14 / melt_fraction
            # tio2 = 3.26 / melt_fraction

            # Kilauea_comparison paper
            # sio2 = 50.19 / melt_fraction
            # al2o3 = 12.79 / melt_fraction
            # feot = 11.34 / melt_fraction
            # mgo = 9.23 / melt_fraction
            # cao = 10.44 / melt_fraction
            # na2o = 2.39 / melt_fraction
            # k2o = 0.43 / melt_fraction
            # p2o5 = 0.27 / melt_fraction
            # mno = 0.18 / melt_fraction
            # tio2 = 2.34 / melt_fraction

            # MORB_comparison paper
            # sio2 = 47.4 / melt_fraction
            # al2o3 = 17.64 / melt_fraction
            # feot = 7.98 / melt_fraction
            # mgo = 7.63 / melt_fraction
            # cao = 12.44 / melt_fraction
            # na2o = 2.65 / melt_fraction
            # k2o = 0.03 / melt_fraction
            # p2o5 = 0.08 / melt_fraction
            # mno = 0.00 / melt_fraction
            # tio2 = 1.01 / melt_fraction

            # Fuego_comparison paper
            sio2 = 51.46 / melt_fraction
            al2o3 = 17.43 / melt_fraction
            feot = 9.42 / melt_fraction
            mgo = 3.78 / melt_fraction
            cao = 7.99 / melt_fraction
            na2o = 3.47 / melt_fraction
            k2o = 0.78 / melt_fraction
            p2o5 = 0.24 / melt_fraction
            mno = 0.19 / melt_fraction
            tio2 = 1.06 / melt_fraction

            # Mt Spurr
            # sio2 = 59.25/ melt_fraction
            # al2o3 = 13.55 / melt_fraction
            # feot = 6.38 / melt_fraction
            # mgo = 3.63 / melt_fraction
            # cao = 4.07 / melt_fraction
            # na2o = 4 / melt_fraction
            # k2o = 1.49 / melt_fraction
            # p2o5 = 0.33 / melt_fraction
            # mno = 0.18 / melt_fraction
            # tio2 = 0.63 / melt_fraction


            

        sio2_n = 100 * sio2 / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        al2o3_n = 100 * al2o3 / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        feot_n = 100 * feot / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        mgo_n = 100 * mgo / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        cao_n = 100 * cao / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        na2o_n = 100 * na2o / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        k2o_n = 100 * k2o / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        p2o5_n = 100 * p2o5 / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        mno_n = 100 * mno / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        tio2_n = 100 * tio2 / (sio2 + al2o3 + feot + mgo + cao + na2o + k2o + p2o5 + mno + tio2)
        self.composition = {"SiO2": sio2_n,
                            "Al2O3": al2o3_n,
                            "TiO2": tio2_n,
                            "FeOT": feot_n,
                            "MgO": mgo_n,
                            "CaO": cao_n,
                            "Na2O": na2o_n,
                            "K2O": k2o_n,
                            "P2O5": p2o5_n,
                            "MnO": mno_n,
                            }

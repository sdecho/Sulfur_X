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

        if choice == 1:
            if melt_fraction > 0.6:
                sio2 = 31.244 * melt_fraction ** 2 - 53.273 * melt_fraction + 72.423
                al2o3 = -2.712 * melt_fraction ** 2 + 4.7856 * melt_fraction + 15.137
                feot = -5.155 * melt_fraction ** 2 + 10.639 * melt_fraction + 5.0263
                mgo = -3.7681 * melt_fraction ** 2 + 8.3599 * melt_fraction + 0.6057
                cao = -6.9226 * melt_fraction ** 2 + 13.182 * melt_fraction + 2.2464
                na2o = 1.0189 * melt_fraction ** 2 - 3.3865 * melt_fraction + 5.2025
                k2o = 0.46459 / melt_fraction
                p2o5 = 0.3592 * melt_fraction ** 2 - 0.6851 * melt_fraction + 0.5388
                mno = 0.0283 * melt_fraction + 0.1444
                tio2 = -0.2832 * melt_fraction + 1.2179
            else:
                melt_fraction = 0.6
                sio2 = 31.244 * melt_fraction ** 2 - 53.273 * melt_fraction + 72.423
                al2o3 = -2.712 * melt_fraction ** 2 + 4.7856 * melt_fraction + 15.137
                feot = -5.155 * melt_fraction ** 2 + 10.639 * melt_fraction + 5.0263
                mgo = -3.7681 * melt_fraction ** 2 + 8.3599 * melt_fraction + 0.6057
                cao = -6.9226 * melt_fraction ** 2 + 13.182 * melt_fraction + 2.2464
                na2o = 1.0189 * melt_fraction ** 2 - 3.3865 * melt_fraction + 5.2025
                k2o = 0.46459 / melt_fraction
                p2o5 = 0.3592 * melt_fraction ** 2 - 0.6851 * melt_fraction + 0.5388
                mno = 0.0283 * melt_fraction + 0.1444
                tio2 = -0.2832 * melt_fraction + 1.2179
        else:
            sio2 = 60.3/ melt_fraction
            al2o3 = 15.8 / melt_fraction
            feot = 5.49 / melt_fraction
            mgo = 2.8 / melt_fraction
            cao = 6.56 / melt_fraction
            na2o = 1.94 / melt_fraction
            k2o = 0.32 / melt_fraction
            p2o5 = 0 / melt_fraction
            mno = 0 / melt_fraction
            tio2 = 0.81 / melt_fraction
            # sio2 = (31.244 - 53.273 + 72.423)/ melt_fraction
            # al2o3 = (-2.712 + 4.7856 + 15.137)/melt_fraction
            # feot = (-5.155 + 10.639 + 5.0263)/melt_fraction
            # mgo = (-3.7681 + 8.3599 + 0.6057)/melt_fraction
            # cao = (-6.9226 + 13.182 + 2.2464)/melt_fraction
            # na2o = (1.0189 - 3.3865 + 5.2025)/melt_fraction
            # k2o = 0.46459 /melt_fraction
            # p2o5 = (0.3592 -0.6851 + 0.5388)/melt_fraction
            # mno = (0.0283 + 0.1444)/melt_fraction
            # tio2 = (-0.2832 + 1.2179)/melt_fraction

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

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
             if melt_fraction > 0.16:
                k2o = 0.5 / melt_fraction
                sio2 = 0.7993 * k2o ** 2 + 1.3594 * k2o + 47.724
                al2o3 = 0.1359 * k2o ** 2 -1.7216 * k2o+ 19.556
                feot = -0.2452 * k2o + 8.3786
                mgo = -0.4446 * k2o ** 2 + 1.1916 * k2o+ 4.2932
                cao = 0.3869 * k2o ** 2 - 4.5494 * k2o + 15.852
                na2o = 0.2359 * k2o ** 2 - 0.2212 * k2o + 3.4682
                p2o5 = -0.0438 * k2o ** 2 + 0.2943 * k2o + 0.5277
                mno = 0.14
                tio2 = 1
             else:
                melt_fraction = 0.16
                k2o = 0.5 / melt_fraction
                sio2 = 0.7993 * k2o ** 2 + 1.3594 * k2o + 47.724
                al2o3 = 0.1359 * k2o ** 2 - 1.7216 * k2o + 19.556
                feot = -0.2452 * k2o + 8.3786
                mgo = -0.4446 * k2o ** 2 + 1.1916 * k2o + 4.2932
                cao = 0.3869 * k2o ** 2 - 4.5494 * melt_fraction + 15.852
                na2o = 0.2359 * k2o ** 2 - 0.2212 * k2o + 3.4682
                p2o5 = -0.0438 * k2o ** 2 + 0.2943 * k2o + 0.5277
                mno = 0.14
                tio2 = 1
            # if melt_fraction > 0.16:
            #     k2o = 0.5 / melt_fraction
            #     sio2 = 4.1544*k2o +45.729
            #     al2o3 = -1.2626*k2o+19.282
            #     feot = 8
            #     mgo = -0.3214*k2o+5.2513
            #     cao = -3.2445*k2o +15.01
            #     na2o = 4
            #     p2o5 = 0.7
            #     tio2 =2
            #     # sio2 = 13.593 * melt_fraction ** 2 - 26.256 * melt_fraction + 56.853
            #     # al2o3 = -11.966 * melt_fraction ** 2 + 16.568 * melt_fraction + 12.54
            #     # feot = -0.9746 * melt_fraction ** 2-0.6548 * melt_fraction + 10.581
            #     # mgo = -3.6759 * melt_fraction ** 2 + 6.8175 * melt_fraction +2.7245
            #     # cao = -16.756 * melt_fraction ** 2 + 24.728 * melt_fraction + 4.08
            #     # na2o = -0.6288 * melt_fraction ** 2 - 0.8552 * melt_fraction + 3.096
            #     # p2o5 = -0.1399* melt_fraction ** 2 - 0.0002 * melt_fraction + 0.197
            #     # mno = -0.5376 * melt_fraction**2 + 0.4872*melt_fraction+0.1881
            #     # tio2 = -1.8215 * melt_fraction**2+1.5844*melt_fraction + 0.8572
            # else:
            #     melt_fraction = 0.16
            #     k2o = 0.5 / melt_fraction
            #     sio2 = 4.1544*k2o +45.729
            #     al2o3 = -1.2626*k2o+19.282
            #     feot = 8
            #     mgo = -0.3214*k2o+5.2513
            #     cao = -3.2445*k2o +15.01
            #     na2o = 4
            #     p2o5 = 0.7
            #     tio2 =2
                # sio2 = 13.593 * melt_fraction ** 2 - 26.256 * melt_fraction + 56.853
                # al2o3 = -11.966 * melt_fraction ** 2 + 16.568 * melt_fraction + 12.54
                # feot = -0.9746 * melt_fraction ** 2-0.6548 * melt_fraction + 10.581
                # mgo = -3.6759 * melt_fraction ** 2 + 6.8175 * melt_fraction +2.7245
                # cao = -16.756 * melt_fraction ** 2 + 24.728 * melt_fraction + 4.08
                # na2o = -0.6288 * melt_fraction ** 2 - 0.8552 * melt_fraction + 3.096
                # k2o =0.24453814354425 / melt_fraction
                # p2o5 = -0.1399* melt_fraction ** 2 - 0.0002 * melt_fraction + 0.197
                # mno = -0.5376 * melt_fraction**2 + 0.4872*melt_fraction+0.1881
                # tio2 = -1.8215 * melt_fraction**2+1.5844*melt_fraction + 0.8572
        else:
            # k2o = 0.5 / melt_fraction
            # sio2 = (0.7993 * 0.5 ** 2 + 1.3594 * 0.5 + 47.724)/melt_fraction
            # al2o3 = (0.1359 * 0.5 ** 2 -1.7216 * 0.5+ 19.556)/melt_fraction
            # feot = (-0.2452 * 0.5 + 8.3786)/melt_fraction
            # mgo = (-0.4446 * 0.5 ** 2 + 1.1916 * 0.5+ 4.2932)/melt_fraction
            # cao = (0.3869 * 0.5 ** 2 - 4.5494 * 0.5 + 15.852)/melt_fraction
            # na2o = (0.2359 * 0.5 ** 2 - 0.2212 * 0.5 + 3.4682)/melt_fraction
            # p2o5 = (-0.0438 * 0.5 ** 2 + 0.2943 * 0.5 + 0.5277)/melt_fraction
            # mno = 0.14/melt_fraction
            # tio2 = 1/melt_fraction
            sio2 = 48.39/ melt_fraction
            al2o3 = 18.83/ melt_fraction
            feot =11.3 / melt_fraction
            mgo = 5.98/ melt_fraction
            cao = 10.10 / melt_fraction
            na2o = 3.21/ melt_fraction
            k2o = 0.43 / melt_fraction
            p2o5 = 0.20/ melt_fraction
            mno = 0.20 / melt_fraction
            tio2 = 1.36/ melt_fraction

            # sio2 = 48.06/ melt_fraction
            # al2o3 = 19.783/ melt_fraction
            # feot =7.686 / melt_fraction
            # mgo = 4.7964/ melt_fraction
            # cao = 13.36375 / melt_fraction
            # na2o = 3.83366/ melt_fraction
            # k2o = 0.8172 / melt_fraction
            # p2o5 = 0.28731/ melt_fraction
            # mno = 0.1375 / melt_fraction
            # tio2 = 1.234/ melt_fraction

            # sio2 = 52.53/ melt_fraction
            # al2o3 = 17.15537/ melt_fraction
            # feot =7.97 / melt_fraction
            # mgo = 4.8283/ melt_fraction
            # cao = 9.5985 / melt_fraction
            # na2o = 3.561047/ melt_fraction
            # k2o = 1.6477 / melt_fraction
            # p2o5 = 0.6628/ melt_fraction
            # mno = 0.127 / melt_fraction
            # tio2 = 1.912/ melt_fraction
            # sio2 = 43.51/ melt_fraction
            # al2o3 = 17.03/ melt_fraction
            # feot =9.026460966 / melt_fraction
            # mgo = 6.023535038/ melt_fraction
            # cao = 12.17536751/ melt_fraction
            # na2o = 1.665036161 / melt_fraction
            # k2o = 0.244538144 / melt_fraction
            # p2o5 = 0.064426649 / melt_fraction
            # mno = 0.1544 / melt_fraction
            # tio2 = 0.6727 / melt_fraction
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

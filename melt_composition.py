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
                pass
            else:
                melt_fraction = 0.16
                
            k2o = 0.5 / melt_fraction
            sio2 = 0.7993 * k2o ** 2 + 1.3594 * k2o + 47.724
            al2o3 = 0.1359 * k2o ** 2 - 1.7216 * k2o + 19.556
            feot = -0.2452 * k2o + 8.3786
            mgo = -0.4446 * k2o ** 2 + 1.1916 * k2o + 4.2932
            cao = 0.3869 * k2o ** 2 - 4.5494 * k2o + 15.852
            na2o = 0.2359 * k2o ** 2 - 0.2212 * k2o + 3.4682
            p2o5 = -0.0438 * k2o ** 2 + 0.2943 * k2o + 0.5277
            mno = 0.14
            tio2 = 3
            
        else: # Input melt composition here if crystallization is disabled.
            sio2 = 51.46/ melt_fraction
            al2o3 = 17.43 / melt_fraction
            feot = 9.42 / melt_fraction
            mgo = 3.78 / melt_fraction
            cao = 7.99 / melt_fraction
            na2o = 3.47 / melt_fraction
            k2o = 0.78 / melt_fraction
            p2o5 = 0.24 / melt_fraction
            mno = 0.19 / melt_fraction
            tio2 = 1.06 / melt_fraction          

        new_composition = {
            "SiO2": sio2,
            "Al2O3": al2o3,
            "TiO2": tio2,
            "FeOT": feot,
            "MgO": mgo,
            "CaO": cao,
            "Na2O": na2o,
            "K2O": k2o,
            "P2O5": p2o5,
            "MnO": mno,
        }
        comp_sum = sum(new_composition.values())
        
        new_composition_normed = {
            k: v*100/comp_sum for k, v in new_composition.items()
        }

        self.composition = new_composition_normed

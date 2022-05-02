import numpy as np
from scipy.optimize import fsolve, root
#  Constants for CO2 solubility on anhydrous base
D_H2O = 2.3
D_ACNK = 3.8
D_FE_MG = -16.3
D_NA_K = 20.1
ALPHA_CO2 = 1
BETA_CO2 = 15.8
C_CO2 = 0.14
B_CO2 = -5.3

# Constants for CO2 solubility on hydrous base
d_H2O = -16.4
d_ACNK = 4.4
d_FE_MG = -17.1
d_NA_K = 22.8
alpha_CO2 = 1
beta_CO2 = 17.3
c_CO2 = 0.12
b_CO2 = -6
# Constants for H2O solubility on anhydrous base
ALPHA_H2O = 0.54
BETA_H2O = 1.24
B_H2O = -2.95
C_H2O = 0.02
# Constants for H2O solubility on hydrous base
alpha_H2O = 0.53
beta_H2O = 2.35
b_H2O = -3.37
c_H2O = -0.02

class IaconoMarziano:
    # This method is coded based on the COH degassing anhydrous-based model from Iacono-Marziano et al. (2012).
    # Please cite the publication above for use of this script
    """ preussure[MPa], temperature_k[K], melt composition: wtsio2[wt% SiO2], wttio2[wt% TiO2], wtal2o3 [wt% al2o3],
    wtfeo[wt% feo,FeO total as FeO], wtmno [wt% MnO], wtmgo[wt% MgO], wtcao[wt% CaO], wtna2o [wt% Na2O],
    wtk2o [wt% k2o], wtp2o5 [wt% p2o5]
    a, b: coefficients for K2O-H2O relation if crystallization is enabled. K2O(wt.%) = a*H2O(wt.%)+b
    """

    def __init__(self, pressure, temperature_k, composition, a, b):
        # Normalize wt % of oxides
        self.Pb = pressure * 10  # transfer from MPa to bar
        self.Tkc = temperature_k

        wtsio2 = composition["SiO2"]
        wttio2 = composition["TiO2"]
        wtal2o3 = composition["Al2O3"]
        wtfeo = composition["FeOT"]
        wtmno = composition["MnO"]
        wtmgo = composition["MgO"]
        wtcao = composition["CaO"]
        wtna2o = composition["Na2O"]
        wtk2o = composition["K2O"]
        wtp2o5 = composition["P2O5"]
        oxide_tot = wtsio2 + wttio2 + wtal2o3 + wtfeo + wtmno + wtmgo + wtcao + wtna2o + wtk2o + wtp2o5
        wtsio2 = wtsio2 / oxide_tot * 100
        wttio2 = wttio2 / oxide_tot * 100
        wtal2o3 = wtal2o3 / oxide_tot * 100
        wtfeo = wtfeo / oxide_tot * 100
        wtmno = wtmno / oxide_tot * 100
        wtmgo = wtmgo / oxide_tot * 100
        wtcao = wtcao / oxide_tot * 100
        wtna2o = wtna2o / oxide_tot * 100
        wtk2o = wtk2o / oxide_tot * 100
        wtp2o5 = wtp2o5 / oxide_tot * 100

        # Convert wt % to mole fractions
        self.nsio2 = wtsio2 / (28.086 + 15.999 * 2)  # moles of siO2
        self.ntio2 = wttio2 / (47.867 + 15.999 * 2)  # moles of tio2
        self.nal2o3 = wtal2o3 / (26.982 * 2 + 15.999 * 3)  # moles of al2o
        self.nfeo = wtfeo / (55.845 + 15.999)  # moles of feo
        self.nmno = wtmno / (54.938 + 15.999)  # moles of mno
        self.nmgo = wtmgo / (24.305 + 15.999)  # moles of mgo
        self.ncao = wtcao / (40.078 + 15.999)  # moles of cao
        self.nna2o = wtna2o / (22.9898 * 2 + 15.999)  # moles of na2o
        self.nk2o = wtk2o / (39.098 * 2 + 15.999)  # moles of k2o
        self.np2o5 = wtp2o5 / (30.973 * 2 + 15.999 * 5)  # moles of p2o5
        # nh2o = wth2o /(15.999 + 2 * 1.0079)  # moles of h2o
        self.ntot = (self.nsio2 + self.ntio2 + self.nal2o3 + self.nfeo + self.nmno + self.nmgo + self.ncao + self.nna2o + self.nk2o + self.np2o5)
        # totalmole
        self.slope_h2o = a  # slope of K2O = a*H2O + b
        self.con_h2o = b  # constant of K2O =a*H2O + b

    def saturation_pressure(self, carbonate_0, h2o_0):
        nh2o = h2o_0 / (15.999 + 2 * 1.0079)
        xh2o = nh2o / (self.ntot + nh2o)
        xsio2 = self.nsio2 / (self.ntot + nh2o)
        xtio2 = self.ntio2 / (self.ntot + nh2o)
        xal2o3 = self.nal2o3 / (self.ntot + nh2o)
        xfeo = self.nfeo / (self.ntot + nh2o)
        xmno = self.nmno / (self.ntot + nh2o)
        xmgo = self.nmgo / (self.ntot + nh2o)
        xcao = self.ncao / (self.ntot + nh2o)
        xna2o = self.nna2o / (self.ntot + nh2o)
        xk2o = self.nk2o / (self.ntot + nh2o)
        xp2o5 = self.np2o5 / (self.ntot + nh2o)
        AI = xal2o3 / (xcao + xna2o + xk2o)
        NBO = 2 * (xh2o + xk2o + xna2o + xcao + xmgo + xfeo - xal2o3) / \
              (2 * xsio2 + 2 * xtio2 + 3 * xal2o3 + xmgo + xfeo + xcao + xna2o + xk2o + xh2o)

        u0 = np.array([self.Pb, 0.5])
        u = root(self.func_initial, u0, (h2o_0, carbonate_0, AI, xfeo + xmgo, xna2o + xk2o, NBO, self.ntot, self.Tkc))
        pressure_sat = u.x[0]
        XH2O_f = u.x[1]
        return pressure_sat, XH2O_f

    def func_initial (self, u, h2o_0, carbonate_0, x_ai, x_feomgo, x_na2ok2o, NBO, anhy_ntot, T):
        P_sat = u[0]
        XH2O_f = u[1]
        eq_H2O = (P_sat ** alpha_H2O) * (XH2O_f ** alpha_H2O) * np.exp(
            beta_H2O * NBO + b_H2O + c_H2O * P_sat / T) - h2o_0
        eq_CO2 = (d_H2O * h2o_0 / (15.999 + 2 * 1.0079) / (anhy_ntot + h2o_0 / (15.999 + 2 * 1.0079))
                  + d_ACNK * x_ai + d_FE_MG * x_feomgo + d_NA_K * x_na2ok2o) \
                 + alpha_CO2 * np.log(P_sat * (1 - XH2O_f)) + beta_CO2 * NBO + b_CO2 + c_CO2 * P_sat / T - np.log(
            carbonate_0)
        F = np.array([eq_H2O, eq_CO2])
        return F

    def coh_solubility(self, Pm, h2o_guess, co2_0, h2o_0, XS_fluid, rS_fluid, u0, choice):

        nh2o = h2o_guess / (15.999 + 2 * 1.0079)
        xh2o = nh2o / (self.ntot + nh2o)
        xsio2 = self.nsio2 / (self.ntot + nh2o)
        xtio2 = self.ntio2 / (self.ntot + nh2o)
        xal2o3 = self.nal2o3 / (self.ntot + nh2o)
        xfeo = self.nfeo / (self.ntot + nh2o)
        xmno = self.nmno / (self.ntot + nh2o)
        xmgo = self.nmgo / (self.ntot + nh2o)
        xcao = self.ncao / (self.ntot + nh2o)
        xna2o = self.nna2o / (self.ntot + nh2o)
        xk2o = self.nk2o / (self.ntot + nh2o)
        xp2o5 = self.np2o5 / (self.ntot + nh2o)
        AI = xal2o3 / (xcao + xna2o + xk2o)
        NBO = 2 * (xh2o + xk2o + xna2o + xcao + xmgo + xfeo - xal2o3) / \
              (2 * xsio2 + 2 * xtio2 + 3 * xal2o3 + xmgo + xfeo + xcao + xna2o + xk2o + xh2o)

        if choice == 1:
            u = root(self.func_crystalization, u0, (h2o_0, co2_0, AI, xfeo + xmgo, xna2o + xk2o, NBO, self.ntot, self.Tkc,
                                 self.Pb, XS_fluid, rS_fluid))

        else:

            u = root(self.func_no_crystallization, u0, (h2o_0, co2_0, AI, xfeo + xmgo, xna2o + xk2o, NBO, self.ntot, self.Tkc,
                               self.Pb, XS_fluid, rS_fluid))
        return u

    def func_crystalization (self, u, h2o_0, co2_0, x_ai, x_feomgo, x_na2ok2o, NBO, anhy_ntot, T, P, XS_fluid, rS_fluid):
        fm = u[0]
        fv = u[1]
        XH2O_f = u[2]
        XCO2_f = u[3]
        H2O_m = u[4]
        CO2_m = u[5]
        fc = u[6]

        eq_H2O = (P ** alpha_H2O) * (XH2O_f ** alpha_H2O) * np.exp(beta_H2O * NBO + b_H2O + c_H2O * P / T) - H2O_m
        eq_CO2 = (d_H2O * (H2O_m / (15.999 + 2 * 1.0079)) / (anhy_ntot + H2O_m / (15.999 + 2 * 1.0079))
                  + d_ACNK * x_ai + d_FE_MG * x_feomgo + d_NA_K * x_na2ok2o) \
                 + alpha_CO2 * np.log(P * XCO2_f) + beta_CO2 * NBO + b_CO2 + c_CO2 * P / T - np.log(CO2_m * 60.009 /44.01)
        eq_totalP = XH2O_f + XCO2_f + XS_fluid - 1
        mb_h2o = fm * H2O_m + fv * 100 * XH2O_f * 18.015 / (XH2O_f * 18.015 + XCO2_f * 44.01 +
                                                            XS_fluid * rS_fluid * 64 + XS_fluid * (
                                                                        1 - rS_fluid) * 34) - h2o_0
        mb_co2 = fm * CO2_m + fv * 1000000 * XCO2_f * 44.01 / (XH2O_f * 18.015 + XCO2_f * 44.01 +
                                                               XS_fluid * rS_fluid * 64 + XS_fluid * (
                                                                           1 - rS_fluid) * 34) - co2_0
        eq_fm_h2o = (h2o_0 * self.slope_h2o + self.con_h2o) / (H2O_m * self.slope_h2o + self.con_h2o) - fm
        eq_mb = fm + fv + fc - 1

        F = np.array([eq_H2O, eq_CO2, eq_totalP, mb_h2o, mb_co2, eq_fm_h2o, eq_mb])
        return F

    def func_no_crystallization (self, u, h2o_0, co2_0, x_ai, x_feomgo, x_na2ok2o, NBO, anhy_ntot, T, P, XS_fluid, rS_fluid):
        fm = u[0]
        fv = u[1]
        XH2O_f = u[2]
        XCO2_f = u[3]
        H2O_m = u[4]
        CO2_m = u[5]

        eq_H2O = (np.log(P * XH2O_f) * alpha_H2O) + (beta_H2O * NBO + b_H2O + c_H2O * P / T) - np.log(H2O_m)
        eq_CO2 = (d_H2O * (H2O_m / (15.999 + 2 * 1.0079)) / (anhy_ntot + H2O_m / (15.999 + 2 * 1.0079))
                  + d_ACNK * x_ai + d_FE_MG * x_feomgo + d_NA_K * x_na2ok2o) \
                 + alpha_CO2 * np.log(P * XCO2_f) + beta_CO2 * NBO + b_CO2 + c_CO2 * P / T - np.log(CO2_m * 60.009 / 44.01)
        eq_totalP = XH2O_f + XCO2_f + XS_fluid - 1
        mb_h2o = fm * H2O_m + fv * 100 * XH2O_f * 18.015 / (XH2O_f * 18.015 + XCO2_f * 44.01 +
                                                            XS_fluid * rS_fluid * 64 + XS_fluid * (
                                                                    1 - rS_fluid) * 34) - h2o_0
        mb_co2 = fm * CO2_m + fv * 1000000 * XCO2_f * 44.01 / (XH2O_f * 18.015 + XCO2_f * 44.01 +
                                                               XS_fluid * rS_fluid * 64 + XS_fluid * (
                                                                       1 - rS_fluid) * 34) - co2_0
        eq_mb = fm + fv - 1

        F = np.array([eq_H2O, eq_CO2, eq_totalP, mb_h2o, mb_co2, eq_mb])
        return F



    def COH_recalc(self, PCO2, PH2O, h2o):
        nh2o = h2o / (15.999 + 2 * 1.0079)
        xh2o = nh2o / (self.ntot + nh2o)
        xsio2 = self.nsio2 / (self.ntot + nh2o)
        xtio2 = self.ntio2 / (self.ntot + nh2o)
        xal2o3 = self.nal2o3 / (self.ntot + nh2o)
        xfeo = self.nfeo / (self.ntot + nh2o)
        xmno = self.nmno / (self.ntot + nh2o)
        xmgo = self.nmgo / (self.ntot + nh2o)
        xcao = self.ncao / (self.ntot + nh2o)
        xna2o = self.nna2o / (self.ntot + nh2o)
        xk2o = self.nk2o / (self.ntot + nh2o)
        xp2o5 = self.np2o5 / (self.ntot + nh2o)
        AI = xal2o3 / (xcao + xna2o + xk2o)
        NBO = 2 * (xh2o + xk2o + xna2o + xcao + xmgo + xfeo - xal2o3) / \
              (2 * xsio2 + 2 * xtio2 + 3 * xal2o3 + xmgo + xfeo + xcao + xna2o + xk2o + xh2o)

        wtH2O = np.exp((np.log(PH2O) * alpha_H2O) + (beta_H2O *NBO + b_H2O + c_H2O * self.Pb / self.Tkc))

        xfeo = self.nfeo/(self.ntot + wtH2O / (15.999 + 2 * 1.0079))
        xmgo = self.nmgo/(self.ntot + wtH2O / (15.999 + 2 * 1.0079))
        xna2o = self.nna2o/(self.ntot + wtH2O / (15.999 + 2 * 1.0079))
        xk2o = self.nk2o/(self.ntot + wtH2O / (15.999 + 2 * 1.0079))
        xal2o3 = self.nal2o3/(self.ntot + wtH2O / (15.999 + 2 * 1.0079))
        xcao = self.ncao/(self.ntot + wtH2O / (15.999 + 2 * 1.0079))
        AI = xal2o3/(xcao + xna2o + xk2o)


        wtcarbonate = np.exp((d_H2O * (wtH2O / (15.999 + 2 * 1.0079)) / (self.ntot + wtH2O / (15.999 + 2 * 1.0079))
                  + d_ACNK * AI + d_FE_MG * (xfeo+xmgo) + d_NA_K * (xna2o + xk2o))\
                 + alpha_CO2 * np.log(PCO2) + beta_CO2 * NBO + b_CO2 + c_CO2 * self.Pb / self.Tkc)
        wtCO2 = 44.01* wtcarbonate/ 60.009
        return wtH2O, wtCO2






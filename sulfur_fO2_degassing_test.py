import numpy as np
from scipy.optimize import fsolve, root
import math

#  Constants for CO2 solubility
D_H2O = 2.3
D_ACNK = 3.8
D_FE_MG = -16.3
D_NA_K = 20.1
ALPHA_CO2 = 1
BETA_CO2 = 15.8
C_CO2 = 0.14
B_CO2 = -5.3

# Constants for H2O solubility on anhydrous base
ALPHA_H2O = 0.54
BETA_H2O = 1.24
B_H2O = -2.95
C_H2O = 0.02

class S_fO2:
    """
    P[MPa],Tkc[K]
    melt composition: wtsio2[wt% SiO2], wttio2[wt% TiO2], wtal2o3 [wt% al2o3],
    wtfeo[wt% feo,FeO total as FeO], wtmno [wt% MnO], wtmgo[wt% MgO], wtcao[wt% CaO], wtna2o [wt% Na2O],
    wtk2o [wt% k2o], wtp2o5 [wt% p2o5], wth2o [wt% h2o]
    model_choice: choice of S-Fe model(0, 100 or other float numbers)
    """

    def __init__(self, pressure, temperature_k, composition, ferric, model_choice):
        # Normalize wt % of oxides
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
        # wth2o = wth2o / oxide_tot * 100

        # Convert wt % to oxide mole fractions
        nsio2 = wtsio2 / (28.086 + 15.999 * 2)  # moles of siO2
        ntio2 = wttio2 / (47.867 + 15.999 * 2)  # moles of tio2
        nal2o3 = wtal2o3 / (26.982 * 2 + 15.999 * 3)  # moles of al2o
        nfeo = wtfeo / (55.845 + 15.999)  # moles of feo
        nmno = wtmno / (54.938 + 15.999)  # moles of mno
        nmgo = wtmgo / (24.305 + 15.999)  # moles of mgo
        ncao = wtcao / (40.078 + 15.999)  # moles of cao
        nna2o = wtna2o / (22.9898 * 2 + 15.999)  # moles of na2o
        nk2o = wtk2o / (39.098 * 2 + 15.999)  # moles of k2o
        np2o5 = wtp2o5 / (30.973 * 2 + 15.999 * 5)  # moles of p2o5
        # nh2o = wth2o /(15.999 + 2 * 1.0079)  # moles of h2o

        # cation mole fraction
        xtot = nsio2 + ntio2+0.5*nal2o3+nfeo+nmno+nmgo+nmgo+0.5*nna2o+0.5*nk2o
        xna = (wtna2o / 30.99) / xtot
        xmg = (wtmgo / 40.32) / xtot
        xal = (wtal2o3 / 50.98) / xtot
        xsi = (wtsio2 / 60.08) / xtot
        xk = (wtk2o / 47.1) / xtot
        xca = (wtcao / 56.08) / xtot
        xti = (wttio2 / 79.9) / xtot
        xmn = (wtmno / 70.94) / xtot
        xfet = (wtfeo / 71.85) / xtot
        xferrous = xfet * (1 - ferric)


        ntot = (nsio2 + ntio2 + nal2o3 + nfeo + nmno + nmgo + ncao + nna2o + nk2o + np2o5)  # totalmole
        self.ntot = ntot
        self.xsio2 = nsio2 / ntot
        self.xtio2 = ntio2 / ntot
        self.xal2o3 = nal2o3 / ntot
        self.xfeo = nfeo / ntot
        self.xmno = nmno / ntot
        self.xmgo = nmgo / ntot
        self.xcao = ncao / ntot
        self.xna2o = nna2o / ntot
        self.xk2o = nk2o / ntot
        self.xp2o5 = np2o5 / ntot
        # self.xh2o = nh2o / ntot
        self.Pb = pressure * 10  # transfer from MPa to bar
        self.Tkc = temperature_k
        self.S_Fe_choice = model_choice # choice of S-Fe model
        # self.NBO = 2 * (self.xk2o + self.xna2o + self.xcao + self.xmgo + self.xfeo -self.xal2o3) / \
        #       (2 * self.xsio2 + 2 * self.xtio2 + 3 * self.xal2o3 + self.xmgo + self.xfeo + self.xcao + self.xna2o + self.xk2o)
        # self.AI = self.xal2o3 / (self.xcao + self.xk2o + self.xna2o)
        # self.h2o_con = ALPHA_H2O * np.log(self.Pb) + BETA_H2O * self.NBO +B_H2O + C_H2O * self.Pb / self.Tkc
        # self.co2_con = D_ACNK * self.AI + D_FE_MG * (self.xfeo + self.xmgo) + D_NA_K * (self.xna2o + self.xk2o) + \
        #                ALPHA_CO2 * np.log(self.Pb) + BETA_CO2 * self.NBO + B_CO2 + C_CO2 * self.Pb / self.Tkc
        self.c_sulfide = 8.77 - 23590 / self.Tkc + (1673 / self.Tkc) * (
                6.7 * (xna + xk) + 4.9 * xmg + 8.1 * xca + 8.9 * (xfet + xmn) + 5 * xti + 1.8 * xal
                - 22.2 * xti * (xfet + xmn) + 7.2 * ((xfet + xmn) * xsi)) - 2.06 * math.erf(-7.2 * (xfet + xmn))
        self.c_sulfate = (-8.02) + (
                21100 + 44000 * xna + 18700 * xmg + 4300 * xal + 35600 * xca + 44200 * xk + 16500 * xferrous + 12600 * xmn) / self.Tkc
        self.lnk = (-55921) / self.Tkc + 25.07 - 0.6465 * np.log(self.Tkc)

    def cohs_solubility(self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons,
                        feo_cr_acc, e_feo_cr, ebalance, u0):
        def func2(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s, feo_cr_acc, e_feo_cr,
                  eb_initial, c_sulfate, lnk, c_sulfide, S_Fe_choice):
            rS_m = u[0]  # S6+/ST in the melt
            rS_f = u[1]  # SO2/ST in the vapor
            rfe3fe2 = u[2]  # Fe3+/Fe2+ in the melt
            fo2 = u[3]  # fO2 in bar

            eq_gas = (fo2 ** 1.5) * ((1 - rS_f) * XS_f * phih2s * P) - (XH2O_f * phih2o * P) * (
                    rS_f * XS_f * phiso2 * P) * (10 ** (4.1245 - 27110 / T))  # P in bar, T in K
            # eq_Fe = np.log(0.5 * rfe3fe2) - (fo2 * 0.196) - fo2_con
            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)

            # relation of S-Fe depending on the choice of S-Fe model
            if S_Fe_choice == 0:
                eq_S_Fe  = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (8.7436*1000000/(T ** 2) - 27703/T + 20.273))
            elif S_Fe_choice == 100:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (- 2863 / T + 7.2))
            elif S_Fe_choice == 1:
                eq_S_Fe = np.exp(c_sulfate-lnk-c_sulfide)*(fo2**2) - (np.exp(c_sulfate-lnk-c_sulfide)*(fo2**2) + 1)*rS_m
            else:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (- 2863 / T + S_Fe_choice))

            eq_e_balance = (wS_m / 10000) * fm * (1 - rS_m) * 8 / 32.065 + \
                           wS_f * fv * ((1 - rS_f) * 8 + 2 * rS_f) / 32.065 + \
                           (1 / (1 + rfe3fe2)) * fm * wFeO_m / (
                                   55.845 + 15.999) + feo_cr_acc + e_feo_cr - eb_initial
            F = np.array([eq_gas, eq_Fe, eq_S_Fe, eq_e_balance])
            return F

        u = root(func2, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                               phi_so2, phi_h2o, phi_h2s,feo_cr_acc, e_feo_cr, ebalance, self.c_sulfate, self.lnk,
                             self.c_sulfide, self.S_Fe_choice), method='lm')
        return u.x

    def cohs_solubility_S6(self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons,feo_cr_acc,
                           e_feo_cr, ebalance, u0):
        def func2(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s, feo_cr_acc, e_feo_cr,
                  eb_initial):

            rS_f = u[0]  # SO2/ST in the vapor
            rfe3fe2 = u[1]  # Fe3+/Fe2+ in the melt
            fo2 = u[2]  # fO2 in bar

            eq_gas = (fo2 ** 1.5) * ((1 - rS_f) * XS_f * phih2s * P) - (XH2O_f * phih2o * P) * (
                    rS_f * XS_f * phiso2 * P) * (10 ** (4.1245 - 27110 / T))  # P in bar, T in K
            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)
            eq_e_balance = wS_f * fv * ((1 - rS_f) * 8 + 2 * rS_f) / 32.065 + \
                           (1 / (1 + rfe3fe2)) * fm * wFeO_m / (55.845 + 15.999) + feo_cr_acc + e_feo_cr - eb_initial
            F = np.array([eq_gas, eq_Fe, eq_e_balance])
            return F

        u1 = root(func2, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                               phi_so2, phi_h2o, phi_h2s, feo_cr_acc, e_feo_cr, ebalance), method='lm')
        return u1.x

    def cohs_solubility_S2(self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons,feo_cr_acc,
                           e_feo_cr, ebalance, u0):
        def func3(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s, feo_cr_acc, e_feo_cr,
                  eb_initial):

            rS_f = u[0]  # SO2/ST in the vapor
            rfe3fe2 = u[1]  # Fe3+/Fe2+ in the melt
            fo2 = u[2]  # fO2 in bar

            eq_gas = (fo2 ** 1.5) * ((1 - rS_f) * XS_f * phih2s * P) - (XH2O_f * phih2o * P) * (
                    rS_f * XS_f * phiso2 * P) * (10 ** (4.1245 - 27110 / T))  # P in bar, T in K
            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)

            eq_e_balance = (wS_m / 10000) * fm * 8 / 32.065 + wS_f * fv * ((1 - rS_f) * 8 + 2 * rS_f) / 32.065 + \
                           (1 / (1 + rfe3fe2)) * fm * wFeO_m / (55.845 + 15.999) + feo_cr_acc + e_feo_cr - eb_initial
            F = np.array([eq_gas, eq_Fe, eq_e_balance])
            return F

        u2 = fsolve(func3, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                                phi_so2, phi_h2o, phi_h2s, feo_cr_acc, e_feo_cr, ebalance))
        return u2

    def cohs_so2(self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons, feo_cr_acc,
                        e_feo_cr, ebalance, u0):
        def func3(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s, feo_cr_acc, e_feo_cr,
                  eb_initial,  c_sulfate, lnk, c_sulfide, S_Fe_choice):
            rS_m = u[0]  # S6+/ST in the melt
            rfe3fe2 = u[1]  # Fe3+/Fe2+ in the melt
            fo2 = u[2]  # fO2 in bar

            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)
            if S_Fe_choice == 0:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (
                            10 ** (8.7436 * 1000000 / (T ** 2) - 27703 / T + 20.273))
            elif S_Fe_choice == 100:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (- 2863 / T + 7.2))
            elif S_Fe_choice == 1:
                eq_S_Fe = np.exp(c_sulfate-lnk-c_sulfide)*(fo2**2) - (np.exp(c_sulfate-lnk-c_sulfide)*(fo2**2) + 1)*rS_m

            else:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (- 2863 / T + S_Fe_choice))
            eq_e_balance = (wS_m / 10000) * fm * (1 - rS_m) * 8 / 32.065 + \
                           wS_f * fv * 2 / 32.065 + \
                           (1 / (1 + rfe3fe2)) * fm * wFeO_m / (
                                   55.845 + 15.999) + feo_cr_acc + e_feo_cr - eb_initial
            F = np.array([eq_Fe, eq_S_Fe, eq_e_balance])
            return F

        u3 = fsolve(func3, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                               phi_so2, phi_h2o, phi_h2s,feo_cr_acc, e_feo_cr, ebalance, self.c_sulfate, self.lnk,
                             self.c_sulfide, self.S_Fe_choice))
        return u3

    def cohs_h2s(self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons,feo_cr_acc,
                        e_feo_cr, ebalance, u0):
        def func4(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s, feo_cr_acc, e_feo_cr,
                  eb_initial, c_sulfate, lnk, c_sulfide, S_Fe_choice):
            rS_m = u[0]  # S6+/ST in the melt
            rfe3fe2 = u[1]  # Fe3+/Fe2+ in the melt
            fo2 = u[2]  # fO2 in bar


            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)
            if S_Fe_choice == 0:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (
                            10 ** (8.7436 * 1000000 / (T ** 2) - 27703 / T + 20.273))
            elif S_Fe_choice == 100:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (- 2863 / T + 7.2))
            elif S_Fe_choice == 1:
                eq_S_Fe = np.exp(c_sulfate-lnk-c_sulfide)*(fo2**2) - (np.exp(c_sulfate-lnk-c_sulfide)*(fo2**2) + 1)*rS_m
            else:
                eq_S_Fe = rS_m * XS_m - (rfe3fe2 ** 8) * (1 - rS_m) * XS_m * (10 ** (- 2863 / T + S_Fe_choice))
            eq_e_balance = (wS_m / 10000) * fm * (1 - rS_m) * 8 / 32.065 + \
                           wS_f * fv * 8/ 32.065 + \
                           (1 / (1 + rfe3fe2)) * fm * wFeO_m / (
                                   55.845 + 15.999) + feo_cr_acc + e_feo_cr - eb_initial
            F = np.array([eq_Fe, eq_S_Fe, eq_e_balance])
            return F

        u4 = fsolve(func4, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                               phi_so2, phi_h2o, phi_h2s, feo_cr_acc, e_feo_cr, ebalance, self.c_sulfate, self.lnk,
                             self.c_sulfide, self.S_Fe_choice))
        return u4

    def cohs_S6_so2(self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons,feo_cr_acc,
                        e_feo_cr, ebalance, u0):
        def func2(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s,feo_cr_acc, e_feo_cr,
                  eb_initial):
            rfe3fe2 = u[0]  # Fe3+/Fe2+ in the melt
            fo2 = u[1]  # fO2 in bar


            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)
            eq_e_balance = wS_f * fv * 2 / 32.065 + (1 / (1 + rfe3fe2)) * fm * wFeO_m / (
                                   55.845 + 15.999) + feo_cr_acc+ e_feo_cr - eb_initial
            F = np.array([eq_Fe, eq_e_balance])
            return F

        u = fsolve(func2, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                               phi_so2, phi_h2o, phi_h2s, feo_cr_acc, e_feo_cr, ebalance))
        return u

    def cohs_S2_h2s (self, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, phi_so2, phi_h2o, phi_h2s, fo2_cons,feo_cr_acc,
                        e_feo_cr, ebalance, u0):
        def func2(u, fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, P, T, fo2_con, phiso2, phih2o, phih2s, feo_cr_acc, e_feo_cr,
                  eb_initial):

            rfe3fe2 = u[0]  # Fe3+/Fe2+ in the melt
            fo2 = u[1]  # fO2 in bar


            eq_Fe = 0.5 * rfe3fe2 - (fo2 ** 0.196) * np.exp(fo2_con)

            eq_e_balance = (wS_m / 10000) * fm * 8 / 32.065 + wS_f * fv * 8 / 32.065 + (1 / (1 + rfe3fe2)) * fm * wFeO_m / (
                                   55.845 + 15.999) + feo_cr_acc+ e_feo_cr - eb_initial
            F = np.array([eq_Fe, eq_e_balance])
            return F

        u = fsolve(func2, u0, (fm, fv, XH2O_f, XS_f, XS_m, wS_m, wS_f, wFeO_m, self.Pb, self.Tkc, fo2_cons,
                               phi_so2, phi_h2o, phi_h2s, feo_cr_acc, e_feo_cr, ebalance))
        return u

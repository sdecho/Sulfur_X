import numpy as np

T0 = 1400  # in Celsius
TK0 = T0 + 273.15  # in kelvin
A = 0.196
B = 11492
C = -6.675
DWFEO = -1.828
DWAL2O3 = -2.243
DWCAO = 3.201
DWNA2O = 5.854
DWK2O = 6.215
E = -3.36
F = -0.000000701  # KPa - 1
G = -1.54e-10  # Pa - 1
H = 3.85e-17  # kPa - 2


class OxygenFugacity:
    """P[MPa],Tkc[K], Initial or average melt composition: wtsio2[wt% SiO2],wttio2[wt% TiO2], wtal2o3 [wt% al2o3],
    wtfeo[wt% feo,FeO total as FeO], wtmno [wt% MnO], wtmgo[wt% MgO], wtcao[wt% CaO], wtna2o [wt% Na2O],
    wtk2o [wt% k2o], wtp2o5 [wt% p2o5]
    logfo2 of fmq buffer, can be calculated with pressure and temperature following Frost(1991)"""

    def __init__(self, P, Tkc, composition):
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

        # Convert wt % to mole fractions
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
        self.ntot = (nsio2 + ntio2 + nal2o3 + nfeo + nmno + nmgo + ncao + nna2o + nk2o + np2o5)  # totalmole

        self.xsio2 = nsio2 / self.ntot
        self.xtio2 = ntio2 / self.ntot
        self.xal2o3 = nal2o3 / self.ntot
        self.xfeo = nfeo / self.ntot
        self.xmno = nmno / self.ntot
        self.xmgo = nmgo / self.ntot
        self.xcao = ncao / self.ntot
        self.xna2o = nna2o / self.ntot
        self.xk2o = nk2o / self.ntot
        self.xp2o5 = np2o5 / self.ntot
        self.Pp = P * 1000000  # transfer from MPa to Pa
        self.Pb = P * 10  # transfer from MPa to bar
        self.Tkc = Tkc
        self.con = B / self.Tkc + C + (
                DWAL2O3 * self.xal2o3 + DWCAO * self.xcao + DWNA2O * self.xna2o + DWK2O * self.xk2o + DWFEO * self.xfeo) \
                   + E * (1 - TK0 / self.Tkc - np.log(self.Tkc / TK0)) + F * self.Pp / self.Tkc + \
                   G * (self.Tkc - TK0) * self.Pp / self.Tkc + H * (self.Pp ** 2) / self.Tkc

    def fe_ratio(self, o2):
        # --------------Kress & Carmichael 1991 parameter values - -------------------
        fo2_ln = np.log(np.power(10, o2))
        rfe_ln = A * fo2_ln + self.con
        rfe = 2 * np.exp(rfe_ln) / (1 + 2 * np.exp(rfe_ln))
        return rfe

    def fo2(self, rfe):
        rfe = (rfe / 2) / (1 - rfe)  # transfer ratio, fe3 / fet to fe3 / fe2;
        o2 = (np.log(rfe) - self.con) / A
        o2 = np.log10(np.exp(o2))
        return o2

    def fmq(self):
        o2_fmq = -25096.3 / self.Tkc + 8.735 + .110 * (self.Pb - 1) / self.Tkc
        return o2_fmq

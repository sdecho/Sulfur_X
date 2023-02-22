import numpy as np
from scipy.stats import truncnorm
R = 8.3145
# T0 = 1400  # in Celsius
# TK0 = T0 + 273.15  # in kelvin
DELV_RXNI = 12.68-22.1
P0_RXNI = 200
T0_RXNI = 1000

P0_RXNIA = 100
T0_RXNIA = 1400
# DELV_RXNII = 45.94-16.764
DELV_RXNII = -49.5
DELH_RXNII = 32379  #delH/R
P0_RXNII = 100  # in MPa
T0_RXNII = 1100 # in C


class PartitionCoefficient:
    """
    P[MPa],Tkc[K], Initial or average melt composition: wtsio2[wt% SiO2], wttio2[wt% TiO2], wtal2o3 [wt% al2o3],
    wtfeo[wt% feo,FeO total as FeO], wtmno [wt% MnO], wtmgo[wt% MgO], wtcao[wt% CaO], wtna2o [wt% Na2O],
    wtk2o [wt% k2o], wtp2o5 [wt% p2o5], wth2o [wt% h2o]
    phih2o, phih2s, and phiso2 : fugacity coefficients of h2o, h2s and so2
    monte (==1) for option to return a random number within the estimated error of each kd
    logfo2 of fmq buffer, can be calculated with pressure and temperature following Frost(1991)
    """

    def __init__(self, P, Tkc, composition, wth2o, phih2o, phih2s, phiso2, monte):
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

        # Normalize wt % of oxides
        oxide_tot = wtsio2 + wttio2 + wtal2o3 + wtfeo + wtmno + wtmgo + wtcao + wtna2o + wtk2o + wtp2o5 + wth2o
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
        wth2o = wth2o / oxide_tot * 100

        # Convert wt % to mole fractions in a single cation base
        nsi = wtsio2 / (28.086 + 15.999 * 2)  # moles of siO2
        nti = wttio2 / (47.867 + 15.999 * 2)  # moles of tio2
        nal = 2 * wtal2o3 / (26.982 * 2 + 15.999 * 3)  # moles of alo1.5
        nfe = wtfeo / (55.845 + 15.999)  # moles of feo
        nmn = wtmno / (54.938 + 15.999)  # moles of mno
        nmg = wtmgo / (24.305 + 15.999)  # moles of mgo
        nca = wtcao / (40.078 + 15.999)  # moles of cao
        nna = 2 * wtna2o / (22.9898 * 2 + 15.999)  # moles of nao0.5
        nk = 2 * wtk2o / (39.098 * 2 + 15.999)  # moles of ko0.5
        nph = 2 * wtp2o5 / (30.973 * 2 + 15.999 * 5)  # moles of po2.5
        self.nh = 2 * wth2o / (15.999 + 2 * 1.0079)  # moles of ho0.5
        self.ntot = (nsi + nti + nal + nfe + nmn + nmg + nca + nna + nk + nph)  # total moles without h2o

        self.xsi = nsi / (self.ntot + self.nh)
        self.xti = nti / (self.ntot + self.nh)
        self.xal = nal / (self.ntot + self.nh)
        self.xfe = nfe / (self.ntot + self.nh)
        self.xmn = nmn / (self.ntot + self.nh)
        self.xmg = nmg / (self.ntot + self.nh)
        self.xca = nca / (self.ntot + self.nh)
        self.xna = nna / (self.ntot + self.nh)
        self.xk = nk / (self.ntot + self.nh)
        self.xph = nph / (self.ntot + self.nh)
        self.xh = self.nh / (self.ntot + self.nh)
        self.Pb = P * 10  # transfer from MPa to bar
        self.Tkc = Tkc
        self.phiso2 = phiso2
        self.ln_feo_coefficient = ((1 - self.xfe) ** 2 * (
                    28870 - 14710 * self.xmg + 1960 * self.xca + 43300 * self.xna + 95380 * self.xk - 76880 * self.xti)
                                   + (1 - self.xfe) * (-62190 * self.xsi + 31520 * (self.xsi ** 2))) / (R * Tkc)
        # self.ln_S_coefficient = 13.406 * self.xh + (TK0/self.Tkc) * (
        #         -8.4995 * self.xfe + 1.8019 * self.xsi - 5.1959 * self.xsi*self.xfe - 10.8588 * self.xca)
        self.ln_S_coefficient = 20.844 * (self.xh ** 2) - 8.7151 * self.xfe  #+ 110.204 * self.xmn
        self.excess_ca = 2 * self.xca - self.xal
        self.excess_na = self.xna - self.xal
        total_oxygen = 2 * (2 * (self.xsi + self.xti) + self.xal * 1.5 + self.xfe + self.xmn + self.xca + 0.5*(self.xna + self.xk + self.xh) + 2.5 * self.xph)
        self.nbo = (total_oxygen - 4 * (self.xsi + self.xti + self.xal + self.xph))/(self.xsi + self.xti + self.xal + self.xph)
        delh_rxn1 = (-0.4739 * (T0_RXNI + 273.15) + 45431) / R
        delh_rxn1a = (-2.7921*(T0_RXNIA + 273.15)-472032) / R

        self.residual_rxn1 = 0.616843 + self.ln_S_coefficient - DELV_RXNI*(self.Pb/10 - P0_RXNI)/(R * self.Tkc) \
                        - (delh_rxn1)*(1/self.Tkc - 1/(T0_RXNI+273.15)) + np.log(self.Pb * phih2o) \
                        - np.log(self.xfe) - self.ln_feo_coefficient - np.log(self.Pb * phih2s)
        # residual_Rxn1 = ln(xh2s) -ln(xS2-) + ln(XH2O_fluid)
        self.residual_rxn1a = 29.79217 + self.ln_S_coefficient - DELV_RXNI*(self.Pb/10 - P0_RXNIA)/(R * self.Tkc) \
                              - (delh_rxn1a) * (1/self.Tkc - 1/(T0_RXNIA+273.15)) - np.log(self.Pb) - np.log(phiso2) \
                              - np.log(self.xfe) - self.ln_feo_coefficient
        # residual_Rxn1a = ln(so2) -ln(xS2-) - 1.5* lnfO2
        # self.residual_rxn2 = -0.2394 - DELV_RXNII*(self.Pb/10 - P0_RXNII)/(R * self.Tkc) \
        #                      - DELH_RXNII*(1/self.Tkc - 1/T0_RXNII) - np.log(phiso2) - np.log(self.Pb)\
        #                      - 8.917 * self.excess_ca - 21.7845 * self.excess_na - 1.7147 * self.nbo
        self.residual_rxn2 = -0.2556 - DELV_RXNII * (self.Pb / 10 - P0_RXNII) / (R * self.Tkc) \
                             - DELH_RXNII * (1 / self.Tkc - 1 / (T0_RXNII+273.15))- np.log(phiso2) - np.log(self.Pb) \
                             - 9.88817 * self.excess_ca - 24.76395 * self.excess_na - 0.97078 * self.nbo
        # self.residual_rxn2 = 0.959178 - DELV_RXNII * (self.Pb / 10 - P0_RXNII) / (R * self.Tkc) \
        #                      - DELH_RXNII * (1 / self.Tkc - 1 / (T0_RXNII + 273.15)) - np.log(phiso2) - np.log(self.Pb) \
        #                      + 2.5883241 * self.excess_ca - 15.0809 * self.excess_na - 0.92322 * self.nbo
        self.monte = monte

    def get_truncated_normal(mean, sd, low, upp):
        return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

    def kd_rxn1(self, xh2o):
        """
        xh2o: mole fraction of h2o in the vapor, 0-1
        """
        rxn1 = self.residual_rxn1 + np.log(xh2o)
        sd_rxn1 = 0.44
        if self.monte == 1:
            rxn1 = np.random.normal(rxn1, sd_rxn1)
        return np.exp(rxn1)

    def kd_rxn1a(self, fo2):
        """
        fo2: oxygen fugacity in bar
        """
        rxn1a = self.residual_rxn1a + 1.5 * np.log(fo2)
        sd_rxn1a = 0.45
        if self.monte == 1:
            rxn1a = np.random.normal(rxn1a, sd_rxn1a)
        return np.exp(rxn1a)

    def kd_rxn2 (self, fo2):
        """
        fo2: oxygen fugacity in bar
        """
        rxn2 = (self.residual_rxn2 - 0.5 * np.log(fo2))
        sd_rxn2 = 0.2136
        if self.monte == 1:
            rxn2 = np.random.normal(rxn2, sd_rxn2)
        return np.exp(rxn2)

    def gas_quilibrium(self, fo2, fh2o, phiso2, phih2s):
        """
        fo2: oxygen fugacity in bar
        fh2o: water fugacity in bar
        phiso2: fugacity coefficient of so2
        phih2s: fugacity coefficient of h2s
        """
        logK = 4.1245 - 27110/self.Tkc
        logfo2 = np.log10(fo2)
        ratio = 10 ** (1.5*logfo2 - logK - np.log10(fh2o)) * phih2s / phiso2
        ratio = ratio / (ratio +1)
        return ratio

    def hydrogen_equilibrium(self, fo2, fh2o):
        logK = 12510/self.Tkc-0.979*np.log10(self.Tkc)+0.483
        logfH2 = np.log10(fh2o)-0.5*np.log10(fo2)-logK
        return 10**logfH2





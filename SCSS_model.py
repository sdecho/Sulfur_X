import numpy as np

# Ideal mixing for SCSS model
B = 9.087
C = -269.40
ASI = -27561
ATI = -11220
AAL = -18450
AMG = -13970
ACA = -7831
AFE = -34274
ANA = -13247
AK = -29015
AH = -17495
ASI_FE = 116568

# P and T for SCSS model


class Sulfur_Saturation:
    """
    P: pressure in MPa
    T: temperature in C
    silicate melt composition: sio2, tio2, al2o3, mgo, mno, feo, cao, na2o, k2o, p2o5, h2o in wt.%
    sulfide composition: Fe, Ni, Cu, O, S in wt.%
    """

    def __init__(self, P, T, composition, h2o, ferric_fe, sulfide_composition):
        self.P = P/1000
        self.Tk = T + 273.15
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
        wth2o = h2o
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

        # Convert wt % to moles on a single cation base
        self.nsio2 = wtsio2 / (28.086 + 15.999 * 2)  # moles of siO2
        self.ntio2 = wttio2 / (47.867 + 15.999 * 2)  # moles of tio2
        self.nal2o3 = 2 * wtal2o3 / (26.982 * 2 + 15.999 * 3)  # moles of alo1/2
        self.nfeot = wtfeo / (55.845 + 15.999)
        self.nfeo = (1 - ferric_fe) * wtfeo / (55.845 + 15.999)  # moles of feo
        self.nfe2o3 = ferric_fe * wtfeo / (55.845 + 15.999)  # moles of feo1/2
        self.nmno = wtmno / (54.938 + 15.999)  # moles of mno
        self.nmgo = wtmgo / (24.305 + 15.999)  # moles of mgo
        self.ncao = wtcao / (40.078 + 15.999)  # moles of cao
        self.nna2o = 2 * wtna2o / (22.9898 * 2 + 15.999)  # moles of na01/2
        self.nk2o = 2 * wtk2o / (39.098 * 2 + 15.999)  # moles of ko1/2
        self.np2o5 = 2 * wtp2o5 / (30.973 * 2 + 15.999 * 5)  # moles of po5/2
        self.nh2o = 2 * wth2o / 18.0153  # moles of ho1/2

        # convert wt.% of sulfide composition to molar fraction
        sulfide_total = sulfide_composition["Fe"] / 55.8457 + sulfide_composition["Ni"] / 58.6934 + \
                        sulfide_composition["Cu"] / 63.546 + sulfide_composition["O"] / 15.9994 + \
                        sulfide_composition["S"] / 32.065
        xfe = 100 *(sulfide_composition["Fe"] / 55.8457) / sulfide_total
        xni = 100*(sulfide_composition["Ni"] / 58.6934) / sulfide_total
        xcu = 100*(sulfide_composition["Cu"]/63.546) / sulfide_total
        xs = 100*(sulfide_composition["S"] / 32.065) / sulfide_total
        xo = 100*(sulfide_composition["O"] / 15.9994) / sulfide_total
        metal_total = xfe +xcu + xni
        # calculate the Fe activity in sulfide
        self.fe_total = xfe / metal_total

    def SCSS_smythe(self):
        # mole fraction of oxides in the silicate melt on a single cation base
        total = self.nsio2 + self.ntio2 + self.nal2o3 + self.nfeo + self.nfe2o3 + self.nmno + self.nmgo + self.ncao + \
                self.nna2o + self.nk2o + self.np2o5 + self.nh2o
        xsio2m = ASI * self.nsio2 / total
        xtio2m = ATI * self.ntio2 / total
        xal2o3m = AAL * self.nal2o3 / total
        xna2om = ANA * self.nna2o / total
        xk2om = AK * self.nk2o / total
        xh2om = AH * self.nh2o / total
        xmgom = AMG * self.nmgo / total
        xfeom = AFE * self.nfeo / total
        xcaom = ACA * self.ncao / total
        xsifem = ASI_FE * (self.nsio2 / total) * (self.nfeo / total)

        # temperature and pressure effect on SCSS
        P_T = self.P/self.Tk
        CP_T = C * P_T
        lnscss = (122175-80.28*self.Tk+8.474*self.Tk*np.log(self.Tk))/(8.314*self.Tk) \
              + B + (xsio2m+xtio2m+xal2o3m+xmgom+xfeom+xcaom+xna2om+xk2om+xh2om+xsifem)/self.Tk+np.log(self.fe_total)-\
              np.log(self.nfeo / total)+CP_T
        scss = np.exp(lnscss) # SCSS in ppm
        return scss

    def SCAS_Zajacz_Tsay(self):
        total_ox = self.nsio2+self.ntio2+0.5*self.nal2o3+self.nfeot+self.nmgo+self.ncao+0.5*self.nna2o+\
                   0.5*self.nk2o+0.5*self.nh2o
        # mole fraction on the oxides base
        xsio2 = self.nsio2/total_ox
        xtio2 = self.ntio2/total_ox
        xal2o3 = 0.5 * self.nal2o3/total_ox
        xfeot = self.nfeot/total_ox
        xmgo = self.nmgo/total_ox
        xcao = self.ncao/total_ox
        xna2o = 0.5 * self.nna2o/total_ox
        xk2o = 0.5 * self.nk2o/total_ox
        xh2o = 0.5 * self.nh2o/total_ox

        #NBO/T
        if (xna2o*2+xk2o*2+2*(xcao+xmgo+xfeot)-xal2o3*2)/(xal2o3*2+xsio2)>0:
            nbot = (xna2o*2+xk2o*2+2*(xcao+xmgo+xfeot)-xal2o3*2)/(xal2o3*2+xsio2)
        else:
            nbot = 0

        #P_Rhyo
        if (xk2o+xcao+xna2o)>xal2o3:
            p_rhyo = 3.11*(xk2o+xcao+xna2o-xal2o3)
        else:
            p_rhyo = 1.54*(xal2o3-(xk2o+xcao+xna2o))

        #P_c
        p_c = ((p_rhyo+251*(xcao**2)+57*(xmgo**2)+154*(xfeot**2))/(2*xal2o3+xsio2))/(1+4.8*nbot)

        #P_T
        p_t = np.exp(-7890/self.Tk)
        #P_H2O
        p_h2o = xh2o*(2.09-1.65*nbot)+0.42*nbot+0.23

        ksp = np.exp(1.226*np.log(p_c*p_t*p_h2o)+0.079)
        xs = ksp/xcao
        scas = xs*total_ox*32.07*10000 # in ppm
        return scas














import numpy as np
import math
from oxygen_fugacity import OxygenFugacity

class Sulfur_Iron:
    def __init__(self, ferric_iron, temperature, model_choice, composition, o2):
        self.ferric = ferric_iron
        self.Tk = temperature + 273.15
        if model_choice == 0:
            self.sulfate = self.Nash()
        elif model_choice == 100:
            self.sulfate = self.Muth()
        elif model_choice == 1:
            self.sulfate = self.OandM(composition, o2)

        else:
            self.sulfate = self.modified(cons=model_choice)

    def Nash(self):
        sulfate_ratio = 10**(8 * np.log10(self.ferric / (1 - self.ferric)) + 8.7436 * 1000000 / (self.Tk**2) -
                             27703 / self.Tk + 20.273)
        sulfate_ratio = sulfate_ratio/(1+sulfate_ratio)
        return sulfate_ratio

    def Muth(self):
        sulfate_ratio =10**(8 * np.log10(self.ferric / (1 - self.ferric)) - 2863 / self.Tk + 7.2)
        sulfate_ratio = sulfate_ratio/(1+sulfate_ratio)
        return sulfate_ratio

    def OandM(self, composition, o2):
        """o2: log10"""
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

        xtot = wtsio2/60.08+wttio2/79.9+wtal2o3/50.98+wtfeo/71.85+wtmgo/40.32+ wtcao/56.08+wtna2o/30.99+wtk2o/47.1+wtmno/70.94

        xna = (wtna2o/30.99)/xtot
        xmg = (wtmgo/40.32)/xtot
        xal = (wtal2o3/50.98)/xtot
        xsi = (wtsio2/60.08)/xtot
        xk = (wtk2o/47.1)/xtot
        xca = (wtcao/56.08)/xtot
        xti = (wttio2/79.9)/xtot
        xmn = (wtmno/70.94)/xtot
        xfet = (wtfeo/71.85)/xtot
        xferrous = xfet * (1 - self.ferric)
        c_sulfide = 8.77-23590/self.Tk+(1673/self.Tk)*(6.7*(xna+xk)+4.9*xmg+8.1*xca+8.9*(xfet+xmn)+5*xti+1.8*xal
                                                       -22.2*xti*(xfet+xmn)+7.2*((xfet+xmn)*xsi))-2.06*math.erf(-7.2*(xfet+xmn))
        c_sulfate = (-8.02) +(21100+44000*xna+18700*xmg+4300*xal+35600*xca+44200*xk+16500*xferrous+12600*xmn)/self.Tk
        lnk = (-55921)/self.Tk+25.07-0.6465*np.log(self.Tk) # SO3/S
        lnrs =(c_sulfate - lnk - c_sulfide) + 2 * np.log(10)*o2
        rs =1-1/(1+np.exp(lnrs))

        return rs

    def modified(self, cons):
        sulfate_ratio = 10 ** (8 * np.log10(self.ferric / (1 - self.ferric)) - 2863 / self.Tk + cons)
        sulfate_ratio = sulfate_ratio / (1 + sulfate_ratio)
        return sulfate_ratio


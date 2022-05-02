import numpy as np


class Sulfur_Iron:
    def __init__(self, ferric_iron, temperature, model_choice):
        self.ferric = ferric_iron
        self.Tk = temperature + 273.15
        if model_choice == 0:
            self.sulfate = self.Nash()
        elif model_choice == 100:
            self.sulfate = self.Muth()
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

    def modified(self, cons):
        sulfate_ratio = 10 ** (8 * np.log10(self.ferric / (1 - self.ferric)) - 2863 / self.Tk + cons)
        sulfate_ratio = sulfate_ratio / (1 + sulfate_ratio)
        return sulfate_ratio


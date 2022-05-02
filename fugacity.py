import numpy as np
P0 = 1  # P0 =1 bar

class Fugacity:
    """Input of Pressure (in MPa) and temperature (in C), output different fugacity coefficients (phiSO2, phiH2S, phiH)"""
    def __init__(self, pressure, temperature):
        self.P = pressure * 10  # pressure in bar
        self.T = temperature + 273.15  # temperature in kelvin
        self.phiSO2 = self.phiso2()
        self.phiH2S = self.phih2s()
        self.phiH2O = self.phih2o()

    def phiso2(self):
        # This function calculates the fugacity coefficient of SO2 following Shi & Saxena 1992. To avoid discontinuity
        # in the result, only the calculation at high pressure is used. Instead, fugacity coefficient is assumed to be 1
        # at pressure lower than 20bar.
        if self.P < 20:
            lnphiSO2 = 0
        else:
            Tcr = 430.95  # critical temperature in K
            Pcr = 78.7295  # critical pressure in bar
            Pr = self.P / Pcr  # reduced pressure
            P0r = P0 / Pcr
            Tr = self.T / Tcr  # reduced temperature

            # Paramatersfor EOS
            AQ1 = 0.92854e00
            AQ2 = 0.43269e-1
            AQ3 = -0.24671e00
            AQ4 = 0
            AQ5 = 0.24999e00
            AQ6 = 0
            AQ7 = -0.53182e00
            AQ8 = -0.16461e-01
            BQ1 = 0.84866e-03
            BQ2 = -0.18379e-02
            BQ3 = 0.66787e-01
            BQ4 = 0
            BQ5 = -0.29427e-01
            BQ6 = 0
            BQ7 = 0.29003e-01
            BQ8 = 0.54808e-02
            CQ1 = -0.35456e-03
            CQ2 = 0.23316e-04
            CQ3 = 0.94159e-03
            CQ4 = 0
            CQ5 = -0.81653e-03
            CQ6 = 0
            CQ7 = 0.23154e-03
            CQ8 = 0.55542e-04

            # ---------EOS from Shi & Saxena 1992 - --------------------------------------
            # Equation (3a) from Shi & Saxena 1992
            A = AQ1 + AQ2 * Tr + AQ3 * (Tr ** -1) + AQ4 * (Tr ** 2) + AQ5 * (Tr ** -2) + AQ6 * (Tr ** 3) + \
                AQ7 * (Tr ** -3) + AQ8 * np.log(Tr)
            B = BQ1 + BQ2 * Tr + BQ3 * (Tr ** -1) + BQ4 * (Tr ** 2) + BQ5 * (Tr ** -2) + BQ6 * (Tr ** 3) + \
                BQ7 * (Tr ** -3) + BQ8 * np.log(Tr)
            C = CQ1 + CQ2 * Tr + CQ3 * (Tr ** -1) + CQ4 * (Tr ** 2) + CQ5 * (Tr ** -2) + CQ6 * (Tr ** 3) + \
                CQ7 * (Tr ** -3) + CQ8 * np.log(Tr)
            # Equation (10) from Shi & Saxena 1992
            Zcomp = A * np.log(Pr / P0r) + B * (Pr - P0r) + (C / 2) * (Pr ** 2 - P0r ** 2)
            # Equation (9) from Shi & Saxena 1992
            lnphiSO2 = Zcomp - np.log(self.P)

        return np.exp(lnphiSO2)

    def phih2s(self):
        # This function calculates the fugacity coefficient of H2S following Shi & Saxena 1992. To avoid discontinuity
        # in the result, only the calculation at high pressure is used. Instead, fugacity coefficient is assumed to be 1
        # at pressure lower than 20bar.
        if self.P < 20:
            lnphiH2S = 0
        else:
            Tcr = 373.55  # critical temperature in K
            Pcr = 90.0779  # critical pressure in bar
            Pr = self.P / Pcr
            P0r = P0 / Pcr
            Tr = self.T / Tcr

            # Table 1 from Shi & Saxena 1992
            AQ1 = 0.59941e00
            AQ2 = -0.15570e-02
            AQ3 = 0.45250e-01
            AQ4 = 0
            AQ5 = 0.36687e00
            AQ6 = 0
            AQ7 = -0.79248e00
            AQ8 = 0.26058e00
            BQ1 = 0.22545e-01
            BQ2 = 0.17473e-02
            BQ3 = 0.48253e-01
            BQ4 = 0
            BQ5 = -0.19890e-01
            BQ6 = 0
            BQ7 = 0.32794e-01
            BQ8 = -0.10985e-01
            CQ1 = 0.57375e-03
            CQ2 = -0.20944e-05
            CQ3 = -0.11894e-02
            CQ4 = 0
            CQ5 = 0.14661e-02
            CQ6 = 0
            CQ7 = -0.75605e-03
            CQ8 = -0.27985e-03

            # Equation (3a) from Shi & Saxena 1992
            A = AQ1 + AQ2 * Tr + AQ3 * (Tr ** -1) + AQ4 * (Tr ** 2) + AQ5 * (Tr ** -2) + AQ6 * (Tr ** 3) + AQ7 * (Tr ** -3) \
                + AQ8 * np.log(Tr)
            B = BQ1 + BQ2 * Tr + BQ3 * (Tr ** -1) + BQ4 * (Tr ** 2) + BQ5 * (Tr ** -2) + BQ6 * (Tr ** 3) + \
                BQ7 * (Tr ** -3) + BQ8 * np.log(Tr)
            C = CQ1 + CQ2 * Tr + CQ3 * (Tr ** -1) + CQ4 * (Tr ** 2) + CQ5 * (Tr ** -2) + CQ6 * (Tr ** 3) + \
                CQ7 * (Tr ** -3) + CQ8 * np.log(Tr)
            # Equation (10) from Shi & Saxena 1992
            Zcomp = A * np.log(Pr / P0r) + B * (Pr - P0r) + (C / 2) * (Pr ** 2 - P0r ** 2)
            # Equation (9) from Shi & Saxena 1992
            lnphiH2S = Zcomp - np.log(self.P)

        return np.exp(lnphiH2S)

    def phih2o(self):
        if self.P < 40:
            lnphi = 0
        else:
            a1 = 2.95177298930e-2
            a2 = -6.33756452413e3
            a3 = -2.75265428882e5
            a4 = 1.29128089283e-3
            a5 = -1.45797416153e2
            a6 = 7.65938947237e4
            a7 = 2.58661493537e-6
            a8 = .52126532146e0
            a9 = -1.39839523753e2
            a10 = -2.36335007175e-8
            a11 = 5.35026383543e-3
            a12 = -.27110649951e0
            a13 = 2.50387836486e4
            a14 = .73226726041e0
            a15 = 1.54833359970e-2
            epsi = 510
            sig = 2.88
            Pm = (3.0636 * sig ** 3 * self.P) / epsi
            Tm = (154 * self.T) / epsi
            R = .08314467 #[cm3 bar / k(mol)]

            # ---------EOS from Zhang & Duan 2009 - --------------------------------------

            Vm = 20  # arbitrary initial guess
            L = 0 # lower bound for Vm
            U = 5 #upper bound for Vm
            P = 0  # just some arbitrary value for code to run
            while abs(Pm - P) > .01:
                Vm = abs(U + L) / 2
                P = (1 + ((a1 + a2 / Tm ** 2 + a3 / Tm ** 3) / Vm) + ((a4 + a5 / Tm ** 2 + a6 / Tm ** 3) / Vm ** 2)
                     + ((a7 + a8 / Tm ** 2 + a9 / Tm ** 3) / Vm ** 4)+ ((a10 + a11 / Tm ** 2 + a12 / Tm ** 3) / Vm ** 5)
                     + ((a13 / (Tm** 3 * Vm ** 2)) * (a14 + (a15 / Vm ** 2)) * np.exp(-a15 / Vm ** 2))) * ((R * Tm) / Vm)
                if P < Pm:
                    U = Vm
                else:
                    L = Vm
            V = 1000 * Vm * (sig / 3.691) ** 3
            # ---------------Equations to find fugacity coefficient - --------------------
            S1 = (a1 + a2 / Tm ** 2 + a3 / Tm ** 3) / Vm + (a4 + a5 / Tm ** 2 + a6 / Tm ** 3) / (2 * Vm ** 2) +\
                 (a7 + a8 / Tm ** 2 + a9 / Tm ** 3) / (4 * Vm ** 4) + (a10 + a11 / Tm ** 2 + a12 / Tm ** 3) / (5 * Vm **5) + \
                 (a13 / (2 * a15 * Tm ** 3)) * (a14 + 1 - (a14 + 1 + (a15 / Vm ** 2)) * np.exp(-a15 / Vm ** 2))
            S2 = (2 * a2 / Tm ** 2 + 3 * a3 / Tm ** 3) / Vm + (2 * a5 / Tm ** 2 + 3 * a6 / Tm ** 3) / 2 * Vm ** 2 + \
                 (2 * a8 / Tm ** 2 + 3 * a9 / Tm ** 3) / 4 * Vm ** 4 + (2 * a11 / Tm ** 2 + 3 * a12 / Tm ** 3) / 5 * Vm ** 5 + \
                 (3 * a13 / 2 * a15 * Tm ** 3) * (a14 + 1 - (a14 + 1 + (a15 / Vm ** 2)) * np.exp(-a15 / Vm ** 2))
            Z = (Pm * Vm) / (R * Tm)
            lnphi = Z - 1 - np.log(Z) + S1
        return np.exp(lnphi)  # fugacity coefficient for pure H2O fluid



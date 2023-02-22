import numpy as np
from scipy.optimize import root
R = 83.14321
MH2O0 = 0.0000328
SIO2 = 50  # SiO2 threshhold for basalt in VolatileCalc is 49wt.%. Here we use 49 wt.% instead.


class VolatileCalc:
    # This script is revised from the published python code MIMic (Rasmussen et al., 2020), which is a direct port of
    # the script developed in the following publication:
    # Newman, Sally, and Jacob B. Lowenstern. "VolatileCalc: a silicate melt–H2O–CO2 solution model written in
    # Visual Basic for excel." Computers & Geosciences 28.5 (2002): 597-604.
    # Please cite both publications above for use of this script

    """ This method is revised from based on the COH degassing anhydrous-based model from Iacono-Marziano et al. (2012).
        preussure[MPa], temperature_k[K], melt composition: wtsio2[wt% SiO2], wttio2[wt% TiO2], wtal2o3 [wt% al2o3],
        wtfeo[wt% feo,FeO total as FeO], wtmno [wt% MnO], wtmgo[wt% MgO], wtcao[wt% CaO], wtna2o [wt% Na2O],
        wtk2o [wt% k2o], wtp2o5 [wt% p2o5]
        a, b: coefficients for K2O-H2O relation if crystallization is enabled. K2O(wt.%) = a*H2O(wt.%)+b
        """

    def __init__(self, TK, sio2, a, b):

        self.Tk = TK
        self.sio2 = sio2
        self.H2Ocoeff = (-0.0000304 + 0.00000129 * self.sio2) / MH2O0
        self.CO2coeff = (0.0000087 - 0.0000001698 * self.sio2) / 0.00000038
        self.slope_h2o = a
        self.constant_h2o =b

    def FNF(self, V,A,B,P):

        return R * self.Tk / (V - B) - A / ((V * V + B * V) * self.Tk**0.5) - P

    def MRK(self, P): #Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
        FNA =(166800000 - 193080 * (self.Tk - 273.15) + 186.4 * (self.Tk - 273.15)**2 - 0.071288 * ((self.Tk - 273.15)**3)) * 1.01325
        FNB = 1.01325 * (73030000 - 71400 * (self.Tk - 273.15) + 21.57 * (self.Tk - 273.15)**2)
        FNC = 1.01325 * (np.exp(
            -11.071 + 5953 / self.Tk - 2746000 / self.Tk ** 2 + 464600000 / self.Tk ** 3) * 0.5 * R * R * self.Tk ** 2.5 / 1.02668 + 40123800)
        B_1 = 14.6
        B_2 = 29.7
        for X_1 in range(2):  # loops twice, once for each CO2 and H2O
            B = X_1 * B_1 + (1 - X_1) * B_2
            A = X_1 ** 2 * FNA + 2 * X_1 * (1 - X_1) * FNC + (1 - X_1) ** 2 * FNB
            Temp2 = B + 5
            Q = 1
            Temp1 = 0
            while abs(Temp2 - Temp1) >= 0.00001:
                Temp1 = Temp2
                F_1 = (self.FNF(Temp1 + 0.01, A, B, P) - self.FNF(Temp1, A, B, P)) / 0.01
                Temp2 = Temp1 - Q * self.FNF(Temp1, A, B, P) / F_1
                F_2 = (self.FNF(Temp2 + 0.01, A, B, P) - self.FNF(Temp2, A, B, P)) / 0.01
                if F_2 * F_1 <= 0:
                    Q = Q / 2.
                if abs(Temp2 - Temp1) > 0.00001:
                    F_1 = F_2
            V = Temp2
            G_1 = np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * FNA + (1 - X_1) * FNC) * np.log(
                (V + B) / V) / (R * self.Tk ** 1.5 * B)
            G_1 = G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * self.Tk ** 1.5 * B ** 2) - np.log(
                P * V / (R * self.Tk))
            G_1 = np.exp(G_1)
            G_2 = np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * FNC + (1 - X_1) * FNB) * np.log(
                (V + B) / V) / (R * self.Tk ** 1.5 * B)
            G_2 = G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * self.Tk ** 1.5 * B ** 2) - np.log(
                P * V / (R * self.Tk))
            G_2 = np.exp(G_2)
            if X_1 == 0:
                fCO2o = G_2 * P #The fugacity of CO2
            if X_1 == 1:
                fH2Oo = G_1 * P #The fugacity of H2O
        return fCO2o, fH2Oo

    def SatPress(self, WtH2O, PPMCO2):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #routine is a list, the first element is M
        #M = 0 (calc saturation pressure), 1 (equilibrium speciation), 2 (isobar calculation)
        #if M = 1 or 2, the second element in M is a starting pressure
        #comp = 0 (basalt), 1 (rhyolite)
        #SiO2 = SiO2 content in weight percent
        #WtH2O = H2O content in weight percent
        #PPMCO2 = CO2 content in ppm
        #TK = temperature in kelvin

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #local variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #If vapor saturation pressure calculation, routine[0] = 1

        Z = 10. #This is the increment of pressure, default is 10 bar increment


        temp = WtH2O + PPMCO2 / 250.
        if temp < 0.5:
            press = 20.
            Z = 1.
        elif temp < 1:
            press = 40.
            Z = 20.
        elif temp < 2:
            press = 100.
            Z = 100.
        elif temp < 3:
            press = 500.
            Z = 100.
        elif temp < 4:
            press = 1000.
            Z = 100.
        elif temp < 5:
            press = 1500.
            Z = 100.
        elif temp < 6:
            press = 2000.
            Z = 100.
        else:
            press = 3000.
            Z = 100.
        changer = press #pressure is the variable
        if self.sio2 <= SIO2:
            comp = 0
        else:
            comp = 1

        #initialize Y for Newton's method
        Y = 1
        GH2O = 2
        GCO2 = 0

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #isobar loop
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        '''
        Uses Newton's method to move changer (WtH2O or P, depending on calling subroutine, and thus M).  Looping stops when
        the two mole fractions add up to 1.
        '''
        [xH2Om, xCO2m, XOHm, xb] = self.water_sp(WtH2O, PPMCO2)

        while abs(GH2O + GCO2 - 1) >= 0.0001:
            #For low-H2O samples, XH2Om is a simple function.  Avoids problems in code resulting in division by zero.
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate xH2Om,xCO2m
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate gas
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Calls MRK procedure to define  the endmember fugacities of H2O and CO2
            fCO2o,fH2Oo = self.MRK(press)
            '''
            GH2O and GCO2 are the mol fractions of H2O and CO2 in the gas phase.  They are dependent on XCO2m and XH2Om
            as well as T and P and endmember H2O and CO2 fugacities.
            Standard state for basalt is 1 bar and 1200C.  All values, plus Delta V and Delta H from Dixon et al. (1995).
            Solubilities at standard state are 0.5 ppm CO2 and 0.11 wt.% H2O for a basalt with 55% SiO2. These values
            are equivalent to mol fractions of .00000038 and .0000328, respectively.
            Solubilities vary accoding to wt.% SiO2 in accord with H2Ocoeff and CO2coeff derived by Dixon (1997) and
            slightly modified so that a basalt with 55% SiO2 is equivalent to the tholeiite in the original (1995) model.

            If GH2O+GCO2 sum to > 1, then pressure must be decreased or either Wt.H2O or PPMCO2 need to be decreased.
            '''

            if comp == 0: #if composition is basaltic
                GH2O = (xH2Om / MH2O0 / self.H2Ocoeff) * (1 / fH2Oo) / np.exp(-12 * (press - 1) / (41.84 * 1.9872 * self.Tk))
                GCO2 = (xCO2m / 0.00000038 / self.CO2coeff) * (1 / fCO2o) / np.exp(-23 * (press - 1) / (41.84 * 1.9872 * self.Tk))
            #
            #Standard state for water in rhyolite is 799 bars and 850C.  Solubility under those conditions is 3.45 wt.% _
            #which corresponds to xH2Om of 0.0323.  These values and partial molar volume of 10 cc/mol from Silver (1977).
            #Delta H of -4420 kcal/mol is for H2O dissolution in albite from Silver et al. (1990).
            #
            #Standard state for CO2 in rhyolite is 750 bars and 850C.  Solubility under those conditions is 340 ppm _
            #which corresponds to xCO2m of 0.000397.  Partial molar volume is 28 cc/mol. Values from Blank et al. (1993).
            #DeltaH for CO2 dissolution in melt (4861 kcal/mol) from Fogel and Rutherford (1990).
            #
            else: #if composition is rhyolitic
                GH2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / np.exp(-5 * (press - 799) / (41.84 * 1.9872 * self.Tk) + 4420 / 1.9872 * (1 / self.Tk - 1 / 1123.16))
                GCO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(-28 * (press - 750) / (41.84 * 1.9872 * self.Tk) + 4861 * (1 / self.Tk - 1 / 1123.16) / 1.9872)

            #Redefine variables
            if abs(GH2O + GCO2 - 1) >= 0.0001:
                Yo = Y
                Y = np.sign(-1. * (GH2O + GCO2 - 1))
                if Y + Yo == 0:
                    Z = Z / 2.

                changer = changer - Z * Y / 2.
                press = changer

        if comp == 0: #if composition is basaltic
            WtH2Om = (xH2Om * 1801.5) / ((xb * 18.015) + (1 - xb) * 36.594) #Calculate Wt% molecular H2O in melt
            WtOHm = (0.5 * XOHm * 1801.5) / ((xb * 18.015) + (1 - xb) * 36.594) #Calculate Wt% hydroxyl in melt
        else: #if composition is rhyolitic
            WtH2Om = (xH2Om * 1801.5) / ((xb * 18.015) + (1 - xb) * 32.5) #Calculate Wt% molecular H2O in melt
            WtOHm = (0.5 * XOHm * 1801.5) / ((xb * 18.015) + (1 - xb) * 32.5) #Calculate Wt% hydroxyl in melt

        if WtH2O < 0.6:
            WtOHm = WtH2O - WtH2Om

        return [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]

    def func2_crys(self, u, P, T, H2Ocoeff, CO2coeff, XS_fluid, rS_fluid, co2_0, h2o_0):

        fm = u[0]
        fv = u[1]
        XH2O_f = u[2]
        XCO2_f = u[3]
        H2O_m = u[4]
        CO2_m = u[5]
        fc = u[6]

        [fCO2o, fH2Oo] = self.MRK(P)
        [xH2Om, xCO2m, XOHm, xb] = self.water_sp(H2O_m, CO2_m)
        # if P/10 >10.7:
        if self.sio2 <= SIO2:  # if composition is basaltic
            eq_H2O = (xH2Om / MH2O0 / H2Ocoeff) * (1 / fH2Oo) / np.exp(
                -12 * (P - 1) / (41.84 * 1.9872 * T)) - XH2O_f
            eq_CO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(
                -23 * (P - 1) / (41.84 * 1.9872 * T)) - XCO2_f
        else:  # if composition is rhyolitic
            eq_H2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / np.exp(
                    -5 * (P - 799) / (41.84 * 1.9872 * T) + 4420 / 1.9872 * (1 / T - 1 / 1123.16)) - XH2O_f
            eq_CO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(
                    -28 * (P - 750) / (41.84 * 1.9872 * T) + 4861 * (1 / T - 1 / 1123.16) / 1.9872) - XCO2_f


            #
            # Standard state for water in rhyolite is 799 bars and 850C.  Solubility under those conditions is 3.45 wt.% _
            # which corresponds to xH2Om of 0.0323.  These values and partial molar volume of 10 cc/mol from Silver (1977).
            # Delta H of -4420 kcal/mol is for H2O dissolution in albite from Silver et al. (1990).
            #
            # Standard state for CO2 in rhyolite is 750 bars and 850C.  Solubility under those conditions is 340 ppm _
            # which corresponds to xCO2m of 0.000397.  Partial molar volume is 28 cc/mol. Values from Blank et al. (1993).
            # DeltaH for CO2 dissolution in melt (4861 kcal/mol) from Fogel and Rutherford (1990).
            #

        # else:
        #     if self.sio2 <= 49:  # if composition is basaltic
        #         eq_CO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(
        #             -23 * (P - 1) / (41.84 * 1.9872 * T)) - XCO2_f
        #     else:  # if composition is rhyolitic
        #         eq_CO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(
        #             -28 * (P - 750) / (41.84 * 1.9872 * T) + 4861 * (1 / T - 1 / 1123.16) / 1.9872) - XCO2_f
        #
        #     # if P / 10 > 10.7:
        #     #     eq_H2O = H2O_m - (0.0003 * P / 10 + 0.2944)
        #     if 2.2 < P / 10 <= 10.7:
        #         eq_H2O = H2O_m - (-0.0007 * (P / 10) **2 + 0.0117 * P / 10 + 0.245)
        #     else:
        #         eq_H2O = H2O_m - (0.0682 * np.log(P / 10) + 0.2032)

        eq_totalP = XH2O_f + XCO2_f + XS_fluid - 1
        mb_h2o = fm * H2O_m + fv * 100 * XH2O_f * 18.015 / (XH2O_f * 18.015 + XCO2_f * 44.01 +
                                                            XS_fluid * rS_fluid * 64 + XS_fluid * (
                                                                        1 - rS_fluid) * 34) - h2o_0
        mb_co2 = fm * CO2_m + fv * 1000000 * XCO2_f * 44.01 / (XH2O_f * 18.015 + XCO2_f * 44.01 +
                                                               XS_fluid * rS_fluid * 64 + XS_fluid * (
                                                                           1 - rS_fluid) * 34) - co2_0
        eq_fm_h2o = (h2o_0 * self.slope_h2o + self.constant_h2o) / (H2O_m * self.slope_h2o + self.constant_h2o) - fm
        eq_mb = fm + fv + fc - 1

        F = np.array([eq_H2O, eq_CO2, eq_totalP, mb_h2o, mb_co2, eq_fm_h2o, eq_mb])
        return F

    def func2_nocrys (self, u, P, T, H2Ocoeff, CO2coeff, XS_fluid, rS_fluid, co2_0, h2o_0):
        fm = u[0]
        fv = u[1]
        XH2O_f = u[2]
        XCO2_f = u[3]
        H2O_m = u[4]
        CO2_m = u[5]

        [fCO2o, fH2Oo] = self.MRK(P)
        [xH2Om, xCO2m, XOHm, xb] = self.water_sp(H2O_m, CO2_m)

        # if P/10 > 10.7:
        if self.sio2 <= SIO2:  # if composition is basaltic
            eq_H2O = (xH2Om / MH2O0 / H2Ocoeff) * (1 / fH2Oo) / np.exp(
                -12 * (P - 1) / (41.84 * 1.9872 * T)) - XH2O_f
            eq_CO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(
                -23 * (P - 1) / (41.84 * 1.9872 * T)) - XCO2_f
        else:  # if composition is rhyolitic
            eq_H2O = (xH2Om/ 0.0323) * (712.2 / fH2Oo) / np.exp(
                    -5 * (P - 799) / (41.84 * 1.9872 * T) + 4420 / 1.9872 * (
                            1 / self.Tk - 1 / 1123.16)) - XH2O_f
            eq_CO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(
                    -28 * (P - 750) / (41.84 * 1.9872 * T) + 4861 * (
                            1 / self.Tk - 1 / 1123.16) / 1.9872) - XCO2_f
            #
            # Standard state for water in rhyolite is 799 bars and 850C.  Solubility under those conditions is 3.45 wt.% _
            # which corresponds to xH2Om of 0.0323.  These values and partial molar volume of 10 cc/mol from Silver (1977).
            # Delta H of -4420 kcal/mol is for H2O dissolution in albite from Silver et al. (1990).
            #
            # Standard state for CO2 in rhyolite is 750 bars and 850C.  Solubility under those conditions is 340 ppm _
            # which corresponds to xCO2m of 0.000397.  Partial molar volume is 28 cc/mol. Values from Blank et al. (1993).
            # DeltaH for CO2 dissolution in melt (4861 kcal/mol) from Fogel and Rutherford (1990).
            #

        # else:
        #     if self.sio2 <= 49:  # if composition is basaltic
        #         eq_CO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(
        #             -23 * (P - 1) / (41.84 * 1.9872 * T)) - XCO2_f
        #     else:  # if composition is rhyolitic
        #         eq_CO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(
        #             -28 * (P - 750) / (41.84 * 1.9872 * T) + 4861 * (1 / T - 1 / 1123.16) / 1.9872) - XCO2_f
        #
        #     # if P / 10 > 10.7:
        #     #     eq_H2O = H2O_m - (0.0003 * P / 10 + 0.2944)
        #     #if 2.2 < P / 10 <= 10.7:
        #         #eq_H2O = H2O_m - (-0.0007 * (P / 10) ** 2 + 0.0117 * P / 10 + 0.245)
        #     #else:
        #     eq_H2O = H2O_m - (0.0682 * np.log(P / 10) + 0.2032)

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

    def coh_solubility(self, Pm, co2_0, h2o_0, XS_fluid, rS_fluid, u0, choice, h2o_guess):
        ''' This method takes input of pressure (MPa), initial H2O (in wt.%), initial CO2 (in ppm),
        sulfur mole fraction in the vapor, initial guess for the solution, and choice of COH degassing model as input.
        This method returns H2O and CO2 contents in the melt and vapor, and mass fractions of vapor, melt and
        crystal (if crystallization is enabled, choice = 1).'''
        if choice == 1:
                u = root(self.func2_crys, u0, (Pm * 10, self.Tk, self.H2Ocoeff, self.CO2coeff, XS_fluid, rS_fluid, co2_0, h2o_0))
        else:
                u = root(self.func2_nocrys, u0, (Pm * 10, self.Tk, self.H2Ocoeff, self.CO2coeff, XS_fluid, rS_fluid, co2_0, h2o_0))


        return u

    def water_sp(self, WtH2O, PPMCO2):

        if WtH2O < 0.6:
            xb = 0
            XOHm = 0
            if self.sio2<=SIO2:  # if the composition is basaltic
                xH2Om = np.exp(-5.827) * WtH2O ** 1.855
                SNCO2 = PPMCO2 * 0.0001 / 44.009  # relative moles of CO2
                SNH2O = WtH2O / 18.015  # relative moles of H2O
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594  # relative moles of oxygen
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)  # mol fraction of CO2 relative to H2O+CO2+O

            else:  # if the composition is rhyolitic
                xH2Om = np.exp(-5.526) * WtH2O ** 1.977
                SNCO2 = PPMCO2 * 0.0001 / 44.009  # relative moles of CO2
                SNH2O = WtH2O / 18.015  # relative moles of H2O
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5  # relative moles of oxygen
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)  # mol fraction of CO2 relative to H2O+CO2+O

        # For WtH2O > 0.6, use method of Silver, Newman, Stolper et al.
        else:  # For WtH2O > 0.6, use method of Silver, Newman, Stolper et al.
            XOHm = 0.01
            XOHmo = XOHm
            SNCO2 = PPMCO2 * 0.0001 / 44.009  # relative moles of CO2
            SNH2O = WtH2O / 18.015  # relative moles of H2O
            if self.sio2<=SIO2:  # if composition is basaltic
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594  # relative moles of oxygen
            else:  # if composition is rhyolitic
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5  # relative moles of oxygen
            XH2O = SNH2O / (SNH2O + SNO)  # mol fraction of water, relative to H2O+O
            XO = 1 - XH2O  # mol fraction of O, relative to H2O+O
            xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)
            if self.sio2<=SIO2:  # if composition is basaltic
                xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 + (
                            PPMCO2 * 0.0001 / 44.009))  # xb is SNH2O over (SNH2O+SNO+SNCO2)
            else:  # if composition is rhyolitic
                xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 + (
                            PPMCO2 * 0.0001 / 44.009))  # xb is SNH2O over (SNH2O+SNO+SNCO2)
            # The following three if/then statements prevent infinite loops and division by zero
            if XH2O - 0.6 * XOHm <= 0:
                XOHm = (0.999 * 2 * XH2O + XOHmo) / 2
            if XO - 0.5 * XOHm <= 0:
                XOHm = (0.999 * 2 * XO + XOHmo) / 2
            if XOHm ** 2 == 0:
                XOHm = XOHmo / 2
            derath = 1
            if self.sio2<=SIO2:  # if composition is basaltic
                while abs(derath) >= 0.00001:  # loop to define XOHm
                    fx = 9.143 - 3.295 * (XOHm - 1) - 2 * 6.019 * (XO - XOHm) - 2 * 0.572 * (XH2O - XOHm) + np.log(
                        XOHm ** 2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                    fxp = -3.295 + 2 * (6.019 + 0.572) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * (
                                (XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (
                                      XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                    derath = fx / fxp
                    XOHm = XOHm - derath
            else:  # if composition is rhyolitic
                while abs(derath) >= 0.00001:  # loop to define XOHm
                    fx = 9.345 - 4.304 * (XOHm - 1) - 2 * 6.277 * (XO - XOHm) - 2 * 2.328 * (XH2O - XOHm) + np.log(
                        XOHm ** 2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                    fxp = -4.304 + 2 * (6.277 + 2.328) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * (
                                (XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (
                                      XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                    derath = fx / fxp
                    XOHm = XOHm - derath
            xH2Om = xb - 0.5 * XOHm
        return xH2Om, xCO2m, XOHm, xb
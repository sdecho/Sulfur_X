import numpy as np
from oxygen_fugacity import OxygenFugacity
from fugacity import Fugacity
from sulfur_partition_coefficients import PartitionCoefficient
from Iacono_Marziano_COH import IaconoMarziano
from melt_composition import MeltComposition
from sulfur_fO2_degassing_test import S_fO2
from VC_COH import VolatileCalc
from newvariables import NewVariables
from S_Fe import Sulfur_Iron
from SCSS_model import Sulfur_Saturation

INC = 20
BAR = 0


class COHS_degassing:
    def __init__(self, pressure, temperature, COH_model, xlt_choice, S_Fe_choice, H2O_initial, CO2_initial,
                 S_initial, a, b, monte_c):
        self.P = pressure
        self.T = temperature
        self.Tk = temperature + 273.15
        self.COH_model = COH_model
        self.xlt_choice = xlt_choice
        self.S_Fe_choice = S_Fe_choice
        self.H2O_0 = H2O_initial
        self.CO2_0 = CO2_initial
        self.S_0 = S_initial
        self.slope_h2o = a
        self.constant_h2o = b
        self.monte = monte_c
        # sulfide composition in wt.%; only relevant if SCSS is of interests.
        self.sulfide = {"Fe": 65.43,
                        "Ni": 0,
                        "Cu": 0,
                        "O": 0,
                        "S": 36.47
                        }

    def definition(self, df_results, index):
        phi_volatiles = Fugacity(self.P, self.T)
        df_results.iloc[index, df_results.columns.get_loc("phi_H2S")] = phi_volatiles.phiH2S
        # df_results["phi_SO2"][1] = phi_volatiles.phiSO2
        # df_results["pressure"][1] = self.P

        return df_results.iloc[index]

    def degassing_redox(self, df_results, index, e_balance_initial, sigma):
        # define the method and calculate the fugacity coefficients of H2O, H2S and SO2
        phi_volatiles = Fugacity(self.P, self.T)
        df_results.iloc[index, df_results.columns.get_loc("phi_H2O")] = phi_volatiles.phiH2O
        df_results.iloc[index, df_results.columns.get_loc("phi_H2S")] = phi_volatiles.phiH2S
        df_results.iloc[index, df_results.columns.get_loc("phi_SO2")] = phi_volatiles.phiSO2
        df_results.iloc[index, df_results.columns.get_loc("pressure")] = self.P

        # define melt composition, COH-only degassing, fo2, and sulfur partition coefficient objects with current
        # pressure i, and melt fraction from previous step
        silicate_melt = MeltComposition(df_results["melt_fraction"][index - 1], self.xlt_choice)

        # define the method calculating the COH degassing
        if self.COH_model == 1:
            coh_degas = VolatileCalc(TK=self.Tk, sio2=silicate_melt.composition["SiO2"],
                                     a=self.slope_h2o, b=self.constant_h2o)

        else:
            coh_degas = IaconoMarziano(pressure=df_results["pressure"][index], temperature_k=self.Tk,
                                       composition=silicate_melt.composition, a=self.slope_h2o, b=self.constant_h2o)

        fo2_degassing = OxygenFugacity(df_results["pressure"][index], self.Tk, silicate_melt.composition)
        re = PartitionCoefficient(df_results["pressure"][index], self.Tk, silicate_melt.composition,
                                  df_results["wH2O_melt"][index - 1],
                                  phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2, self.monte)
        Fe2_cr = 0
        Fe3_cr = 0

        # absolute log10fO2 is recalculated with a fixed delta_FMQ
        # df_results["fO2"][index] = fo2_degassing.fmq() + delta_FMQ
        # Fe3+/FeT is recalculated with the fixed delta_FMQ
        rs_melt = Sulfur_Iron(ferric_iron=df_results["ferric_ratio"][index - 1], temperature=self.T,
                              model_choice=self.S_Fe_choice, composition=silicate_melt.composition,
                              o2=df_results["fO2"][index - 1])

        def_variables = NewVariables(df_results["pressure"][index],
                                     300)  # the inputs of P and l are not essential here.
        fo2_tr, XH2O_fluid_tr, XCO2_fluid_tr, XSO2_f_tr, XH2S_f_tr, XS_f_tr, ferric_ratio_tr, wS_f_tr, XS6_m_tr, \
        XS2_m_tr, XS_m_tr, wH2O_m_tr, wCO2_m_tr, wS_m_tr, kd1_tr, kd2_tr, kd1a_tr, kd_combined_tr, rs_m_tr, \
        rs_f_tr, melt_fraction_tr, crystal_fraction_tr, vapor_fraction_tr \
            = def_variables.iteration_v(df_results["XH2O_fluid"][index - 1], df_results["ferric_ratio"][index - 1],
                                        df_results["wH2O_melt"][index - 1], df_results["wCO2_melt"][index - 1],
                                        df_results["wS_melt"][index - 1])
        # if fO2 tracker is enabled, fO2 change due to S degassing and crystallization of Fe3+ and Fe2+
        # initial fO2 of current step using Fe3+/FeT ratio from the previous step
        fo2_tr[0] = 10 ** (fo2_degassing.fo2(df_results["ferric_ratio"][index - 1]))
        # water fugacity using the XH2O from previous P step
        fH2O_initial = df_results["XH2O_fluid"][index - 1] * df_results["pressure"][index] * phi_volatiles.phiH2O * 10
        rs_fluid_initial = re.gas_quilibrium(fo2=fo2_tr[0], fh2o=fH2O_initial, phiso2=phi_volatiles.phiSO2,
                                             phih2s=phi_volatiles.phiH2S)
        # initial XSO2/XST in the fluid using fO2 and fH2O from previous step
        rs_f_tr[0] = rs_fluid_initial
        XH2O_fluid_tr[0] = df_results["XH2O_fluid"][index - 1]
        XCO2_fluid_tr[0] = df_results["XCO2_fluid"][index - 1]
        XSO2_f_tr[0] = df_results["XSO2_fluid"][index - 1]
        XH2S_f_tr[0] = df_results["XH2S_fluid"][index - 1]
        wH2O_m_tr[0] = df_results["wH2O_melt"][index - 1]
        wCO2_m_tr[0] = df_results["wCO2_melt"][index - 1]
        wS_m_tr[0] = df_results["wS_melt"][index - 1]
        kd1_tr[0] = re.kd_rxn1(xh2o=XH2O_fluid_tr[0])
        kd2_tr[0] = re.kd_rxn2(fo2=fo2_tr[0])
        kd1a_tr[0] = re.kd_rxn1a(fo2=fo2_tr[0])
        ferric_ratio_tr[0] = df_results["ferric_ratio"][index - 1]
        rs_melt_initial = rs_melt.sulfate
        rs_m_tr[0] = rs_melt_initial
        XS_m_tr[0] = (df_results["wS_melt"][index - 1] / (10000 * 32.065)) / \
                     (re.ntot + df_results["wS_melt"][index - 1] / (10000 * 32.065) + re.nh +
                      df_results["wCO2_melt"][index - 1] / (10000 * 44.01))
        XS6_m_tr[0] = rs_melt_initial * XS_m_tr[0]
        XS2_m_tr[0] = (1 - rs_melt_initial) * XS_m_tr[0]

        if self.P >= BAR:
            # kd_combined_tr[0] = (rs_fluid_initial * kd1a_tr[0] + (1 - rs_fluid_initial) * kd1_tr[0]) * (
            #         1 - rs_melt_initial) + rs_melt_initial * kd2_tr[0]
            kd_combined_tr[0] = kd1_tr[0] * (1 - rs_melt_initial) + rs_melt_initial * kd2_tr[0]
        else:
            kd_combined_tr[0] = df_results["kd_combined_molar"][index - 1] + INC

        XS_f_tr[0] = XS_m_tr[0] * kd_combined_tr[0]
        melt_fraction_tr[0] = df_results["melt_fraction"][index - 1]
        vapor_fraction_tr[0] = df_results["vapor_fraction"][index - 1]
        crystal_fraction_tr[0] = df_results["crystal_fraction"][index - 1]

        for j in range(1, def_variables.n):
            # calculate initial solution for the system from the previous iteration
            XH2S_f_initial = kd_combined_tr[j - 1] * (1 - rs_f_tr[j - 1]) * XS_m_tr[j - 1]
            XSO2_f_initial = kd_combined_tr[j - 1] * rs_f_tr[j - 1] * XS_m_tr[j - 1]
            XH2O_f_initial = XH2O_fluid_tr[j - 1]
            XCO2_f_initial = 1 - XSO2_f_initial - XH2S_f_initial - XH2O_f_initial

            # recalculate H2O, CO2 in the melt, melt fraction, vapor fraction and crystal fraction after S degassing
            if self.xlt_choice == 1:  # crystallization is enabled
                initial_guess = np.array(
                    [melt_fraction_tr[j - 1], vapor_fraction_tr[j - 1], XH2O_f_initial,
                     XCO2_f_initial, wH2O_m_tr[j - 1], wCO2_m_tr[j - 1], crystal_fraction_tr[j - 1]])
                root = coh_degas.coh_solubility(Pm=self.P, h2o_guess=df_results["wH2O_melt"][index - 1],
                                                co2_0=self.CO2_0,
                                                h2o_0=self.H2O_0, XS_fluid=XSO2_f_initial + XH2S_f_initial,
                                                rS_fluid=rs_f_tr[j - 1], u0=initial_guess, choice=self.xlt_choice)

                if root.x[6] <= 0:
                    crystal_fraction_tr[j] = 0
                else:
                    crystal_fraction_tr[j] = root.x[6]
                if root.x[1] <= 0:
                    vapor_fraction_tr[j] = 0
                else:
                    vapor_fraction_tr[j] = root.x[1]

            else:
                initial_guess = np.array(
                    [melt_fraction_tr[j - 1], vapor_fraction_tr[j - 1], XH2O_f_initial, XCO2_f_initial,
                     wH2O_m_tr[j - 1], wCO2_m_tr[j - 1]])
                root = coh_degas.coh_solubility(Pm=self.P, h2o_guess=df_results["wH2O_melt"][index - 1],
                                                co2_0=self.CO2_0,
                                                h2o_0=self.H2O_0, XS_fluid=XSO2_f_initial + XH2S_f_initial,
                                                rS_fluid=rs_f_tr[j - 1], u0=initial_guess, choice=self.xlt_choice)
                crystal_fraction_tr[j] = 0
                if root.x[1] <= 0:
                    vapor_fraction_tr[j] = 0
                else:
                    vapor_fraction_tr[j] = root.x[1]

            fm = 1 - vapor_fraction_tr[j] - crystal_fraction_tr[j]

            fv = vapor_fraction_tr[j]
            XH2O_f = root.x[2]
            XCO2_f = root.x[3]
            wtH2O_m = root.x[4]
            wtCO2_m = root.x[5]
            # update melt composition, fm, fv (and fc), H2O, CO2 in the melt with new values
            melt_comp_updated = MeltComposition(melt_fraction=fm, choice=self.xlt_choice)
            melt_fraction_tr[j] = fm
            XH2O_fluid_tr[j] = XH2O_f
            XCO2_fluid_tr[j] = XCO2_f
            wH2O_m_tr[j] = wtH2O_m
            wCO2_m_tr[j] = wtCO2_m
            fo2_degassing = OxygenFugacity(self.P, self.Tk, melt_comp_updated.composition)
            re_update = PartitionCoefficient(self.P, self.Tk, melt_comp_updated.composition, wtH2O_m,
                                             phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2,
                                             self.monte)

            if self.xlt_choice == 1:  # redox budget of the crystals assuming crystals always take the same Fe3+/FeT as the melt in the previous step
                Fe3_cr = df_results["ferric_ratio"][index - 1] * (
                        df_results["FeOT"][index - 1] * df_results["melt_fraction"][index - 1] - fm *
                        melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
                Fe2_cr = (1 - df_results["ferric_ratio"][index - 1]) * \
                         (df_results["FeOT"][index - 1] * df_results["melt_fraction"][index - 1] - fm *
                          melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
            else:  # if crystallization is disabled
                Fe2_cr = 0
                Fe3_cr = 0

            wS_f_tr[j] = 100 * (XSO2_f_initial + XH2S_f_initial) * 32.065 / (
                    XH2O_f * 18.015 + XCO2_f * 44.01 + XSO2_f_initial * 64 + XH2S_f_initial * 34)
            kdS_wt = wS_f_tr[j] * 10000 / df_results["wS_melt"][index - 1]

            if vapor_fraction_tr[j] + crystal_fraction_tr[j] > 0:
                DS_wt = kdS_wt * vapor_fraction_tr[j] / (vapor_fraction_tr[j] + crystal_fraction_tr[j])
            else:
                DS_wt = kdS_wt
            wtS_m = self.S_0 / (fm * (1 - DS_wt) + DS_wt)
            XS_m_tr[j] = (wtS_m / (10000 * 32.065)) / (
                    re_update.ntot + wtS_m / (10000 * 32.065) + 2 * wtH2O_m / 18.015 + wtCO2_m / (10000 * 44.01))
            XS6_m_tr[j] = XS_m_tr[j] * rs_m_tr[j - 1]
            XS2_m_tr[j] = XS_m_tr[j] * (1 - rs_m_tr[j - 1])

            cohs = S_fO2(self.P, self.Tk, melt_comp_updated.composition, ferric=ferric_ratio_tr[j-1], model_choice=self.S_Fe_choice)
            ## calculate the fO2 at the current iteration using Nash (Fe3+/FeT, S6+/ST), gas equilibrium (fO2, SO2/ST),
            # Kress and Carmicheal (Fe3+/FeT, fO2) and redox budget conservation
            ## Initial guess of S6+/ST, SO2/ST, Fe3+/FeT and fO2 are from previous iteration.
            low_l = 0.0001
            up_l = 0.9999

            if low_l < rs_m_tr[j - 1] < up_l:
                if low_l < rs_f_tr[j - 1] < up_l:

                    initial_guess_2 = np.array(
                        [rs_m_tr[j - 1], rs_f_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]),
                         fo2_tr[j - 1]])

                    [rs_m_n, rs_f_n, ferric_ratio_n, fo2_n] = \
                        cohs.cohs_solubility(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                             XS_m=XS_m_tr[j], wS_f=wS_f_tr[j], wS_m=wS_m_tr[j],
                                             wFeO_m=melt_comp_updated.composition["FeOT"],
                                             phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                             phi_h2s=phi_volatiles.phiH2S,
                                             fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                             e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                             u0=initial_guess_2)

                elif rs_f_tr[j - 1] > up_l:
                    initial_guess_2 = np.array(
                        [rs_m_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                    [rs_m_n, ferric_ratio_n, fo2_n] = \
                        cohs.cohs_so2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                      XS_m=XS_m_tr[j],
                                      wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                      phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                      phi_h2s=phi_volatiles.phiH2S,
                                      fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                      e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                      u0=initial_guess_2)
                    rs_f_n = 1
                else:
                    initial_guess_2 = np.array(
                        [rs_m_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]),
                         fo2_tr[j - 1]])
                    [rs_m_n, ferric_ratio_n, fo2_n] = \
                        cohs.cohs_h2s(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                      XS_m=XS_m_tr[j],
                                      wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                      phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                      phi_h2s=phi_volatiles.phiH2S,
                                      fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                      e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                      u0=initial_guess_2)
                    rs_f_n = 0

            elif rs_m_tr[j - 1] > up_l:
                if rs_f_tr[j - 1] > up_l:
                    initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                    [ferric_ratio_n, fo2_n] = \
                        cohs.cohs_S6_so2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                         XS_m=XS_m_tr[j],
                                         wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                         phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                         phi_h2s=phi_volatiles.phiH2S,
                                         fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                         e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                         u0=initial_guess_2)
                    rs_f_n = 1
                    rs_m_n = 1
                elif low_l < rs_f_tr[j - 1] < up_l:
                    initial_guess_2 = np.array(
                        [rs_f_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])
                    [rs_f_n, ferric_ratio_n, fo2_n] = \
                        cohs.cohs_solubility_S6(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                                XS_m=XS_m_tr[j],
                                                wS_f=wS_f_tr[j], wS_m=wS_m_tr[j],
                                                wFeO_m=melt_comp_updated.composition["FeOT"],
                                                phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                                phi_h2s=phi_volatiles.phiH2S,
                                                fo2_cons=fo2_degassing.con,
                                                feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                                e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                                u0=initial_guess_2)
                    rs_m_n = 1
                else:
                    initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                    [ferric_ratio_n, fo2_n] = \
                        cohs.cohs_S6_so2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                         XS_m=XS_m_tr[j],
                                         wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                         phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                         phi_h2s=phi_volatiles.phiH2S,
                                         fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                         e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                         u0=initial_guess_2)
                    rs_f_n = 0
                    rs_m_n = 1
            else:
                if rs_f_tr[j - 1] < low_l:
                    initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                    [ferric_ratio_n, fo2_n] = \
                        cohs.cohs_S2_h2s(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                         XS_m=XS_m_tr[j],
                                         wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                         phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                         phi_h2s=phi_volatiles.phiH2S,
                                         fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                         e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                         u0=initial_guess_2)
                    rs_f_n = 0
                    rs_m_n = 0
                elif low_l < rs_f_tr[j - 1] < up_l:
                    initial_guess_2 = np.array(
                        [rs_f_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])
                    [rs_f_n, ferric_ratio_n, fo2_n] = \
                        cohs.cohs_solubility_S2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                                XS_m=XS_m_tr[j],
                                                wS_f=wS_f_tr[j], wS_m=wS_m_tr[j],
                                                wFeO_m=melt_comp_updated.composition["FeOT"],
                                                phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                                phi_h2s=phi_volatiles.phiH2S,
                                                fo2_cons=fo2_degassing.con,
                                                feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                                e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                                u0=initial_guess_2)
                    rs_m_n = 0
                else:
                    initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                    [ferric_ratio_n, fo2_n] = \
                        cohs.cohs_S2_h2s(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                         XS_m=XS_m_tr[j],
                                         wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                         phi_so2=phi_volatiles.phiSO2, phi_h2o=phi_volatiles.phiH2O,
                                         phi_h2s=phi_volatiles.phiH2S,
                                         fo2_cons=fo2_degassing.con, feo_cr_acc=df_results["ferrous_cr"][index - 1],
                                         e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                         u0=initial_guess_2)
                    rs_f_n = 1
                    rs_m_n = 0

            ferric_ratio_tr[j] = ferric_ratio_n / (ferric_ratio_n + 1)
            fo2_tr[j] = fo2_n
            rs_f_tr[j] = rs_f_n
            rs_m_tr[j] = rs_m_n

            kd1_tr[j] = re_update.kd_rxn1(xh2o=XH2O_fluid_tr[j])
            kd2_tr[j] = re_update.kd_rxn2(fo2=fo2_tr[j])
            kd1a_tr[j] = re_update.kd_rxn1a(fo2=fo2_tr[j])
            if self.P >= BAR:
                kd_combined_n = (kd1_tr[j] * (1 - rs_f_n) + kd1a_tr[j] * rs_f_n) * (1 - rs_m_n) + kd2_tr[j] * rs_m_n
            else:
                kd_combined_n = df_results["kd_combined_molar"][index - 1] + INC
            kd_combined_tr[j] = kd_combined_n
            XS_f_tr[j] = kd_combined_tr[j] * df_results["XS_melt"][index - 1]
            XSO2_f_tr[j] = XS_f_tr[j] * rs_f_tr[j]
            XH2S_f_tr[j] = XS_f_tr[j] * (1 - rs_f_tr[j])

            wS_f_tr[j] = 100 * (XSO2_f_tr[j] + XH2S_f_tr[j]) * 32.065 / (
                    XH2O_fluid_tr[j] * 18.015 + XCO2_fluid_tr[j] * 44.01 + XSO2_f_tr[j] * 64 + XH2S_f_tr[j] * 34)

            kdS_wt = wS_f_tr[j] * 10000 / df_results["wS_melt"][index - 1]

            if vapor_fraction_tr[j] + crystal_fraction_tr[j] > 0:
                DS_wt = kdS_wt * vapor_fraction_tr[j] / (vapor_fraction_tr[j] + crystal_fraction_tr[j])
            else:
                DS_wt = kdS_wt
            wtS_m = self.S_0 / (melt_fraction_tr[j] * (1 - DS_wt) + DS_wt)
            wS_m_tr[j] = wtS_m
            XS_m_tr[j] = (wtS_m / (10000 * 32.065)) / (
                    re_update.ntot + wtS_m / (10000 * 32.065) + 2 * wtH2O_m / 18.015 + wtCO2_m / (10000 * 44.01))
            XS6_m_tr[j] = XS_m_tr[j] * rs_m_tr[j]
            XS2_m_tr[j] = XS_m_tr[j] * (1 - rs_m_tr[j])
            if j > 30:
                print(j)

            # if abs(XS_m_tr[j]-XS_m_tr[j-1]) < 0.0001:
            if abs(np.log10(fo2_tr[j]) - np.log10(fo2_tr[j - 1])) < sigma:
                df_results.iloc[index, df_results.columns.get_loc("kd_RxnI")] = kd1_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_RxnIa")] = kd1a_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_RxnII")] = kd2_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("SO2/ST")] = rs_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_combined_molar")] = kd_combined_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_combined_wt")] = kdS_wt
                df_results.iloc[index, df_results.columns.get_loc("DS_bulk")] = DS_wt
                df_results.iloc[index, df_results.columns.get_loc("XS_melt")] = XS_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("wS_melt")] = wS_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("wS_fluid")] = wS_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XS_fluid")] = XH2S_f_tr[j] + XSO2_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XH2S_fluid")] = XH2S_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XSO2_fluid")] = XSO2_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("SO2_fugacity")] = XSO2_f_tr[j] * phi_volatiles.phiSO2 * self.P * 10
                df_results.iloc[index, df_results.columns.get_loc("H2S_fugacity")] = XH2S_f_tr[j] * phi_volatiles.phiH2S * self.P * 10
                df_results.iloc[index, df_results.columns.get_loc("XH2O_fluid")] = XH2O_fluid_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XCO2_fluid")] = XCO2_fluid_tr[j]
                if wH2O_m_tr[j] <= df_results["wH2O_melt"][index - 1]:
                    df_results.iloc[index, df_results.columns.get_loc("wH2O_melt")] = wH2O_m_tr[j]
                else:
                    df_results.iloc[index, df_results.columns.get_loc("wH2O_melt")] = df_results["wH2O_melt"][index - 1]
                df_results.iloc[index, df_results.columns.get_loc("wCO2_melt")] = wCO2_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("melt_fraction")] = melt_fraction_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("vapor_fraction")] = vapor_fraction_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("crystal_fraction")] = crystal_fraction_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("water_fugacity")] = df_results["XH2O_fluid"][
                                                          index] * self.P * 10 * phi_volatiles.phiH2O
                df_results.iloc[index, df_results.columns.get_loc("fO2")] = np.log10(fo2_tr[j])
                df_results.iloc[index, df_results.columns.get_loc("FMQ")] = fo2_degassing.fmq()
                df_results.iloc[index, df_results.columns.get_loc("ferric_ratio")] = ferric_ratio_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("S6+/ST")] = rs_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("FeOT")] = melt_comp_updated.composition["FeOT"]
                # df_results["Fe_cr"][i] = df_results["Fe_cr"][i-1]+e_FeO_cr

                break
            elif j == def_variables.n:
                df_results.iloc[index, df_results.columns.get_loc("kd_RxnI")] = kd1_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_RxnIa")] = kd1a_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_RxnII")] = kd2_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("SO2/ST")] = rs_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_combined_molar")] = kd_combined_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("kd_combined_wt")] = kdS_wt
                df_results.iloc[index, df_results.columns.get_loc("DS_bulk")] = DS_wt
                df_results.iloc[index, df_results.columns.get_loc("XS_melt")] = XS_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("wS_melt")] = wS_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("wS_fluid")] = wS_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XS_fluid")] = XH2S_f_tr[j] + XSO2_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XH2S_fluid")] = XH2S_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XSO2_fluid")] = XSO2_f_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("SO2_fugacity")] = XSO2_f_tr[j] * phi_volatiles.phiSO2 * self.P * 10
                df_results.iloc[index, df_results.columns.get_loc("H2S_fugacity")] = XH2S_f_tr[j] * phi_volatiles.phiH2S * self.P * 10
                df_results.iloc[index, df_results.columns.get_loc("XH2O_fluid")] = XH2O_fluid_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("XCO2_fluid")] = XCO2_fluid_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("wH2O_melt")] = wH2O_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("wCO2_melt")] = wCO2_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("melt_fraction")] = melt_fraction_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("vapor_fraction")] = vapor_fraction_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("crystal_fraction")] = crystal_fraction_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("water_fugacity")] = df_results["XH2O_fluid"][
                                                          index] * self.P * 10 * phi_volatiles.phiH2O
                df_results.iloc[index, df_results.columns.get_loc("fO2")] = fo2_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("ferric_ratio")] = ferric_ratio_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("S6+/ST")] = rs_m_tr[j]
                df_results.iloc[index, df_results.columns.get_loc("FeOT")] = melt_comp_updated.composition["FeOT"]
                # df_results["Fe_cr"][i] = df_results["Fe_cr"][i-1]+e_FeO_cr
                df_results.iloc[index, df_results.columns.get_loc("FMQ")] = fo2_degassing.fmq()
                print("Calculation do not converge")
        melt_comp_updated = MeltComposition(df_results["melt_fraction"][index], choice=self.xlt_choice)
        re_new = PartitionCoefficient(self.P, self.Tk, melt_comp_updated.composition, df_results["wH2O_melt"][index],
                                      phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2,
                                      self.monte)
        solubility = Sulfur_Saturation(P=self.P, T=self.T, composition=melt_comp_updated.composition,
                                       h2o=df_results["wH2O_melt"][index], ferric_fe=df_results["ferric_ratio"][index],
                                       sulfide_composition=self.sulfide)
        df_results.iloc[index, df_results.columns.get_loc("sulfate_m")] = df_results["melt_fraction"][index] * (df_results["wS_melt"][index] / 10000) * \
                                         df_results["S6+/ST"][index] / 32.065
        df_results.iloc[index, df_results.columns.get_loc("sulfide_m")] = df_results["melt_fraction"][index] * (df_results["wS_melt"][index] / 10000) * \
                                         (1 - df_results["S6+/ST"][index]) / 32.065
        df_results.iloc[index, df_results.columns.get_loc("SO2_f")] = df_results["vapor_fraction"][index] * df_results["wS_fluid"][index] * \
                                     df_results["SO2/ST"][index] / 32.065
        df_results.iloc[index, df_results.columns.get_loc("H2S_f")] = df_results["vapor_fraction"][index] * df_results["wS_fluid"][index] * (
                1 - df_results["SO2/ST"][index]) / 32.065
        df_results.iloc[index, df_results.columns.get_loc("ferric")] = melt_comp_updated.composition["FeOT"] * df_results["melt_fraction"][index] * \
                                      df_results["ferric_ratio"][index] / (55.845 + 15.999)
        df_results.iloc[index, df_results.columns.get_loc("ferrous")] = melt_comp_updated.composition["FeOT"] * df_results["melt_fraction"][index] * \
                                       (1 - df_results["ferric_ratio"][index]) / (55.845 + 15.999)
        df_results.iloc[index, df_results.columns.get_loc("ferric_cr")] = df_results["ferric_cr"][index - 1] + Fe3_cr
        df_results.iloc[index, df_results.columns.get_loc("ferrous_cr")] = df_results["ferrous_cr"][index - 1] + Fe2_cr
        df_results.iloc[index, df_results.columns.get_loc("electron_balance")] = df_results["sulfide_m"][index] * 8 + df_results["H2S_f"][index] * 8 \
                                                + 2 * df_results["SO2_f"][index] + df_results["ferrous"][index] + \
                                                df_results["ferrous_cr"][index]
        df_results.iloc[index, df_results.columns.get_loc("SCSS")] = solubility.SCSS_smythe()
        df_results.iloc[index, df_results.columns.get_loc("SCAS")] = solubility.SCAS_Zajacz_Tsay()
        df_results.iloc[index, df_results.columns.get_loc("fH2")] = re_new.hydrogen_equilibrium(fo2=10 ** df_results["fO2"][index],
                                                               fh2o=df_results["water_fugacity"][index])

        return df_results.iloc[index]

    def degassing_noredox(self, delta_FMQ, df_results, index):
        # redefine the fugacity coefficient object at new pressure i, only related to P and T
        phi_volatiles = Fugacity(self.P, self.T)
        df_results.iloc[index, df_results.columns.get_loc("phi_H2O")] = phi_volatiles.phiH2O
        df_results.iloc[index, df_results.columns.get_loc("phi_H2S")] = phi_volatiles.phiH2S
        df_results.iloc[index, df_results.columns.get_loc("phi_SO2")] = phi_volatiles.phiSO2
        df_results.iloc[index, df_results.columns.get_loc("pressure")] = self.P
        # define melt composition, COH-only degassing, fo2, and sulfur partition coefficient objects with current pressure i, and melt fraction from previous step
        silicate_melt = MeltComposition(df_results["melt_fraction"][index - 1], self.xlt_choice)

        if self.COH_model == 1:
            coh_degas = VolatileCalc(TK=self.Tk, sio2=silicate_melt.composition["SiO2"], a=self.slope_h2o,
                                     b=self.constant_h2o)

        else:
            coh_degas = IaconoMarziano(pressure=df_results["pressure"][index], temperature_k=self.Tk,
                                       composition=silicate_melt.composition, a=self.slope_h2o, b=self.constant_h2o)

        fo2_degassing = OxygenFugacity(df_results["pressure"][index], self.Tk, silicate_melt.composition)
        re = PartitionCoefficient(df_results["pressure"][index], self.Tk, silicate_melt.composition,
                                  df_results["wH2O_melt"][index - 1],
                                  phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2, self.monte)

        # absolute log10fO2 is recalculated with a fixed delta_FMQ
        df_results.iloc[index, df_results.columns.get_loc("fO2")] = fo2_degassing.fmq() + delta_FMQ
        # Fe3+/FeT is recalculated with the fixed delta_FMQ
        df_results.iloc[index, df_results.columns.get_loc("ferric_ratio")] = fo2_degassing.fe_ratio(df_results["fO2"][index])
        rs_melt = Sulfur_Iron(ferric_iron=df_results["ferric_ratio"][index], temperature=self.T,
                              model_choice=self.S_Fe_choice, composition= silicate_melt.composition, o2=df_results["fO2"][index - 1])

        ## calculate three kds, SO2/ST in the vapor using the water fugacity from previous step, and fO2 from this pressure step.
        df_results.iloc[index, df_results.columns.get_loc("kd_RxnI")] = re.kd_rxn1(xh2o=df_results["XH2O_fluid"][index - 1])
        df_results.iloc[index, df_results.columns.get_loc("kd_RxnII")] = re.kd_rxn2(fo2=10 ** (df_results["fO2"][index]))
        df_results.iloc[index, df_results.columns.get_loc("kd_RxnIa")] = re.kd_rxn1a(fo2=10 ** df_results["fO2"][index])
        df_results.iloc[index, df_results.columns.get_loc("S6+/ST")] = rs_melt.sulfate
        df_results.iloc[index, df_results.columns.get_loc("SO2/ST")] = re.gas_quilibrium(fo2=10 ** df_results["fO2"][index],
                                                        fh2o=df_results["water_fugacity"][index - 1],
                                                        phiso2=phi_volatiles.phiSO2, phih2s=phi_volatiles.phiH2S)
        # combined molar Kd weighed by SO2/ST in the vapor and S6+/ST in the melt
        df_results.iloc[index, df_results.columns.get_loc("kd_combined_molar")] = df_results["kd_RxnI"][index] * \
                                                                                  (1 - df_results["S6+/ST"][index]) + \
                                                                                  df_results["S6+/ST"][index] * df_results["kd_RxnII"][index]
        ## S mole fraction in the vapor using combined molar Kd and mole fraction of S in the melt from previous step
        df_results.iloc[index, df_results.columns.get_loc("XS_fluid")] = df_results["XS_melt"][index - 1] * df_results["kd_combined_molar"][index]
        df_results.iloc[index, df_results.columns.get_loc("XSO2_fluid")] = df_results["XS_fluid"][index] * df_results["SO2/ST"][index]
        df_results.iloc[index, df_results.columns.get_loc("XH2S_fluid")] = df_results["XS_fluid"][index] * (1 - df_results["SO2/ST"][index])

        if self.xlt_choice == 1:  # if crystallization is enabled
            # With known S contents in the melt, using COH model and mass balance to solve for the CO2, H2O in the melt
            # and in the vapor, and mass fractions of vapor, melt and crystal
            # initial guess comes from the previous degassing step
            initial_guess = np.array(
                [df_results["melt_fraction"][index - 1], df_results["vapor_fraction"][index - 1],
                 df_results["XH2O_fluid"][index - 1], df_results["XCO2_fluid"][index - 1],
                 df_results["wH2O_melt"][index - 1],
                 df_results["wCO2_melt"][index - 1], df_results["crystal_fraction"][index - 1]])
            root = coh_degas.coh_solubility(Pm=df_results["pressure"][index],
                                            h2o_guess=df_results["wH2O_melt"][index - 1], co2_0=self.CO2_0,
                                            h2o_0=self.H2O_0, XS_fluid=df_results["XS_fluid"][index],
                                            rS_fluid=df_results["SO2/ST"][index],
                                            u0=initial_guess, choice=self.xlt_choice)

            if root.x[1] <= 0:
                df_results.iloc[index, df_results.columns.get_loc("vapor_fraction")] = 0
            else:
                df_results.iloc[index, df_results.columns.get_loc("vapor_fraction")] = root.x[1]

            if root.x[6] <= 0:
                df_results.iloc[index, df_results.columns.get_loc("crystal_fraction")] = 0
            else:
                df_results.iloc[index, df_results.columns.get_loc("crystal_fraction")] = root.x[6]

        else:  # if crystallization is disabled
            initial_guess = np.array(
                [df_results["melt_fraction"][index - 1], df_results["vapor_fraction"][index - 1],
                 df_results["XH2O_fluid"][index - 1], df_results["XCO2_fluid"][index - 1],
                 df_results["wH2O_melt"][index - 1],
                 df_results["wCO2_melt"][index - 1]])
            root = coh_degas.coh_solubility(Pm=df_results["pressure"][index],
                                            h2o_guess=df_results["wH2O_melt"][index - 1],
                                            co2_0=self.CO2_0,
                                            h2o_0=self.H2O_0, XS_fluid=df_results["XS_fluid"][index],
                                            rS_fluid=df_results["SO2/ST"][index],
                                            u0=initial_guess, choice=self.xlt_choice)
            if root.x[1] <= 0:
                df_results.iloc[index, df_results.columns.get_loc("vapor_fraction")] = 0

            else:
                df_results.iloc[index, df_results.columns.get_loc("vapor_fraction")] = root.x[1]

            df_results.iloc[index, df_results.columns.get_loc("crystal_fraction")] = 0

        fm = 1 - df_results["vapor_fraction"][index] - df_results["crystal_fraction"][index]

        # fv = root.x[1]
        XH2O_f = root.x[2]
        XCO2_f = root.x[3]
        wtH2O_m = root.x[4]
        wtCO2_m = root.x[5]

        ## update melt composition after calculating the melt fraction fm
        melt_comp_updated = MeltComposition(melt_fraction=fm, choice=self.xlt_choice)
        ## update mass fractions of melt, vapor, CO2, H2O contents in the vapor and melt
        df_results.iloc[index, df_results.columns.get_loc("melt_fraction")] = fm
        df_results.iloc[index, df_results.columns.get_loc("XH2O_fluid")] = XH2O_f
        df_results.iloc[index, df_results.columns.get_loc("XCO2_fluid")] = XCO2_f
        df_results.iloc[index, df_results.columns.get_loc("wH2O_melt")] = wtH2O_m
        # if wtH2O_m <= df_results["wH2O_melt"][i-1]:
        #     df_results["wH2O_melt"][i] = wtH2O_m
        # else:
        #     df_results["wH2O_melt"][i] = df_results["wH2O_melt"][i-1]
        df_results.iloc[index, df_results.columns.get_loc("wCO2_melt")] = wtCO2_m

        fo2_degassing = OxygenFugacity(df_results["pressure"][index], self.Tk, melt_comp_updated.composition)
        re_update = PartitionCoefficient(df_results["pressure"][index], self.Tk, melt_comp_updated.composition, wtH2O_m,
                                         phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2, self.monte)
        df_results.iloc[index, df_results.columns.get_loc("wS_fluid")] = 100 * (
                df_results["XSO2_fluid"][index] + df_results["XH2S_fluid"][index]) * 32.065 / (
                                                XH2O_f * 18.015 + XCO2_f * 44.01 + df_results["XSO2_fluid"][
                                            index] * 64 + df_results["XH2S_fluid"][index] * 34)
        df_results.iloc[index, df_results.columns.get_loc("kd_combined_wt")] = df_results["wS_fluid"][index] * 10000 / df_results["wS_melt"][index - 1]

        if df_results["vapor_fraction"][index] + df_results["crystal_fraction"][index] > 0:
            df_results.iloc[index, df_results.columns.get_loc("DS_bulk")] = df_results["kd_combined_wt"][index] * df_results["vapor_fraction"][index] / (
                    df_results["vapor_fraction"][index] + df_results["crystal_fraction"][index])
        else:
            df_results.iloc[index, df_results.columns.get_loc("DS_bulk")] = df_results["kd_combined_wt"][index]

        df_results.iloc[index, df_results.columns.get_loc("wS_melt")] = self.S_0 / (
                fm * (1 - df_results["DS_bulk"][index]) + df_results["DS_bulk"][index])
        df_results.iloc[index, df_results.columns.get_loc("XS_melt")] = (df_results["wS_melt"][index] / (10000 * 32.065)) / \
                                       (re_update.ntot + df_results["wS_melt"][index] / (10000 * 32.065) +
                                        df_results["wH2O_melt"][index] / 18.015 +
                                        df_results["wCO2_melt"][index] / (10000 * 44.01))
        df_results.iloc[index, df_results.columns.get_loc("water_fugacity")] = df_results["XH2O_fluid"][index] * df_results["pressure"][index] \
                                              * phi_volatiles.phiH2O * 10
        df_results.iloc[index, df_results.columns.get_loc("SO2_fugacity")] = df_results["XSO2_fluid"][index] * phi_volatiles.phiSO2 * \
                                            df_results["pressure"][index] * 10
        df_results.iloc[index, df_results.columns.get_loc("H2S_fugacity")] = df_results["XH2S_fluid"][index] * phi_volatiles.phiH2S * \
                                            df_results["pressure"][index] * 10

        if self.xlt_choice == 1:  ## redox budget of the crystals assuming crystals always take the same Fe3+/FeT as the melt in the previous step
            Fe3_cr = df_results["ferric_ratio"][index] * (
                    df_results["FeOT"][index] * df_results["melt_fraction"][index] -
                    fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
            Fe2_cr = (1 - df_results["ferric_ratio"][index]) * (
                    df_results["FeOT"][index] * df_results["melt_fraction"][index]
                    - fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
            # e_FeO_cr = (1 - df_results["ferric_ratio"][i - 1]) * \
            #            (df_results["FeOT"][i-1] * df_results["melt_fraction"][i-1] - fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
        else:  ## if crystallization is disabled, =0
            Fe2_cr = 0
            Fe3_cr = 0

        melt_comp_updated = MeltComposition(df_results["melt_fraction"][index], choice=self.xlt_choice)
        re_new = PartitionCoefficient(self.P, self.Tk, melt_comp_updated.composition, df_results["wH2O_melt"][index],
                                      phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2,
                                      self.monte)
        solubility = Sulfur_Saturation(P=self.P, T=self.T, composition=melt_comp_updated.composition,
                                       h2o=df_results["wH2O_melt"][index], ferric_fe=df_results["ferric_ratio"][index],
                                       sulfide_composition=self.sulfide)
        df_results.iloc[index, df_results.columns.get_loc("sulfate_m")] = df_results["melt_fraction"][index] * (df_results["wS_melt"][index] / 10000) * \
                                         df_results["S6+/ST"][index] / 32.065
        df_results.iloc[index, df_results.columns.get_loc("sulfide_m")] = df_results["melt_fraction"][index] * (df_results["wS_melt"][index] / 10000) * \
                                         (1 - df_results["S6+/ST"][index]) / 32.065
        df_results.iloc[index, df_results.columns.get_loc("SO2_f")] = df_results["vapor_fraction"][index] * df_results["wS_fluid"][index] * \
                                     df_results["SO2/ST"][index] / 32.065
        df_results.iloc[index, df_results.columns.get_loc("H2S_f")] = df_results["vapor_fraction"][index] * df_results["wS_fluid"][index] * (
                1 - df_results["SO2/ST"][index]) / 32.065
        df_results.iloc[index, df_results.columns.get_loc("ferric")] = melt_comp_updated.composition["FeOT"] * df_results["melt_fraction"][index] * \
                                      df_results["ferric_ratio"][index] / (55.845 + 15.999)
        df_results.iloc[index, df_results.columns.get_loc("ferrous")] = melt_comp_updated.composition["FeOT"] * df_results["melt_fraction"][index] * \
                                       (1 - df_results["ferric_ratio"][index]) / (55.845 + 15.999)
        df_results.iloc[index, df_results.columns.get_loc("ferric_cr")] = df_results["ferric_cr"][index - 1] + Fe3_cr
        df_results.iloc[index, df_results.columns.get_loc("ferrous_cr")] = df_results["ferrous_cr"][index - 1] + Fe2_cr
        df_results.iloc[index, df_results.columns.get_loc("electron_balance")] = df_results["sulfide_m"][index] * 8 + df_results["H2S_f"][index] * 8 \
                                                + 2 * df_results["SO2_f"][index] + df_results["ferrous"][index] + \
                                                df_results["ferrous_cr"][index]
        df_results.iloc[index, df_results.columns.get_loc("SCSS")] = solubility.SCSS_smythe()
        df_results.iloc[index, df_results.columns.get_loc("SCAS")] = solubility.SCAS_Zajacz_Tsay()
        df_results.iloc[index, df_results.columns.get_loc("fH2")] = re_new.hydrogen_equilibrium(
            fo2=10 ** df_results["fO2"][index],
            fh2o=df_results["water_fugacity"][index])

        return df_results.iloc[index]

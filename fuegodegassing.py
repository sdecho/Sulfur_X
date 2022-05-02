import pandas as pd
import numpy as np
from oxygen_fugacity import OxygenFugacity
from fugacity import Fugacity
from sulfur_partition_coefficients import PartitionCoefficient
import matplotlib.pyplot as plt
from Iacono_Marziano_COH import IaconoMarziano
from melt_composition import MeltComposition
from sulfur_fO2_degassing_test import S_fO2
from VC_COH import VolatileCalc
from newvariables import NewVariables

df = pd.read_csv("Fuego.csv")
df_2 = pd.read_csv("Hawaii.csv")
## Input of all the initial conditions
temperature = 1030 #float(input("What is the temperature in celsius?")) # temperature in C
tk = float(temperature) + 273.15  # in kelvin
ferric_ratio_0 = 0.215 #float(input("What is the initial Ferric/total Fe ratio in the melt?")) # Fe3+/FeT
H2O_initial = 4.5  #float(input("What is the initial water contents in the melt in %?")) ## wt.%
CO2_initial = 3300 #float(input("What is the initial carbon dioxide contents in the melt in ppm?")) # CO2 initial in ppm
carbonate_initial = CO2_initial * 60.009 / 44.01
S_initial = 2650 #float(input("What is the initial sulfur contents in the melt in ppm?"))  # in ppm
choice = 1 #int(input("Crystallization(=1) or no crystallization(0)?")) # choice ==1, crystallization
COH_model =0 #int(input("Which COH degassing model to use? 1 for VolatileCalc; 0 for Iacono_Marziano"))  # if =1, VolcatileCalc; =0, Iacono_marziano model
fo2_tracker =1 #int(input("fO2 change with S degassing? =1, yes; =0, no."))  # if =0, do not change fO2; = 1, fo2 changes

# Initiate composition module to calculate melt composition. If crystallization is enabled, melt composition is a function of degree of melting. If not, melt composition is constant
melt_comp_initial = MeltComposition(melt_fraction=1, choice=choice)
# Calculating the initial vapor saturation pressure, XH2O, XCO2 using initial water and CO2.
m_run = 10
df_S_m = pd.DataFrame()
df_CO2_m = pd.DataFrame()

for i in range (0, m_run):


    if COH_model == 1: ## if VolatileCalc is chosen
        VC_coh = VolatileCalc(tk, melt_comp_initial.composition["SiO2"])
        [P_initial, WtH2O, PPMCO2, WtH2Om, WtOHm, XH2Of_initial, XCO2] = VC_coh.SatPress(WtH2O=H2O_initial,
                                                                                         PPMCO2=CO2_initial)
        print(f" The initial vapor saturation pressure is {P_initial} bar, and the initial vapor concentration is XH2O = {XH2Of_initial}"
            f" and XCO2 = {1 - XH2Of_initial}")
    else:## if IaconoMarziano model is chosen
        coh = IaconoMarziano(600, tk, melt_comp_initial.composition)
        [P_initial, XH2Of_initial] = coh.saturation_pressure(carbonate_initial, H2O_initial)  # P_initial in bar
        print(
            f" The initial vapor saturation pressure is {P_initial} bar, and the initial vapor concentration is XH2O = {XH2Of_initial}"
            f" and XCO2 = {1 - XH2Of_initial}.")
        coh = IaconoMarziano((P_initial) / 10, tk, melt_comp_initial.composition)
    # Define the variables
    l = 300 # The total steps of decompression degassing
    def_variables = NewVariables(P_initial, l)
    my_data = def_variables.results_dic()
    df_results = pd.DataFrame(data=my_data)


    # # Define vectors to store values through fO2-SCOH degassing iteration
    fo2_tr, XH2O_fluid_tr, XCO2_fluid_tr, XSO2_f_tr, XH2S_f_tr, XS_f_tr, ferric_ratio_tr, wS_f_tr, XS6_m_tr, \
    XS2_m_tr, XS_m_tr, wH2O_m_tr, wCO2_m_tr, wS_m_tr, kd1_tr, kd2_tr, kd1a_tr, kd_combined_tr, rs_m_tr, \
    rs_f_tr, melt_fraction_tr, crystal_fraction_tr, vapor_fraction_tr \
        = def_variables.iteration_v(XH2Of_initial, ferric_ratio_0, H2O_initial, CO2_initial, S_initial)

    # # initiate OxygenFugacity module to calculate fO2 with given ferric iron ratio or vice versa
    fo2_0 = OxygenFugacity(P_initial / 10, tk, melt_comp_initial.composition)
    # initiate Fugacity module to calculate fugaicyt coefficients for H2O, SO2, H2S
    phi = Fugacity(P_initial / 10, temperature)

    re = PartitionCoefficient(P_initial / 10, tk, melt_comp_initial.composition, H2O_initial, phi.phiH2O, phi.phiH2S,
                              phi.phiSO2)

    ## Initial value of different variables
    XS_initial = (S_initial / (10000 * 32.065)) / (
            re.ntot + S_initial / (10000 * 32.065) + re.nh + CO2_initial / (10000 * 44.01)) # Initial sulfur mole fraction in the melt
    fH2O_initial = XH2Of_initial * P_initial * phi.phiH2O ## initial water fugacity calculated by initial water contents in the vapor
    rs_melt_initial = 10**(8 * np.log10(ferric_ratio_0 / (1 - ferric_ratio_0)) - 2863/tk + 6.8)
    # rs_melt_initial = 10**(8 * np.log10(ferric_ratio_0 / (1 - ferric_ratio_0)) + 8.7436 * 1000000/(tk ** 2) - 27703/tk + 20.273)
                              ## initial S6+/S2- in the melt calculated by Nash et al .(2019) model
    rs_melt_initial = rs_melt_initial / (1 + rs_melt_initial) ## initial S6+/ST
    e_balance_initial = (S_initial / 10000) * (1 - rs_melt_initial) * 8 / 32.065 \
                        + (1 - ferric_ratio_0) * melt_comp_initial.composition["FeOT"] / (55.845 + 15.999) ## moles of RB (redox budget) relative to S6+ and Fe3+ in 100 g of starting material
    df_results["wS_melt"][0] = S_initial # in ppm
    df_results["wH2O_melt"][0] = H2O_initial # in wt.%
    df_results["wCO2_melt"][0] = CO2_initial ## in ppm
    df_results["XS_melt"][0] = XS_initial ## mole fraction in the melt
    df_results["fO2"][0] = fo2_0.fo2(ferric_ratio_0)  ## log10fO2
    df_results["XCO2_fluid"][0] = 1 - XH2Of_initial ## mole fraction, in the vapor
    df_results["XH2O_fluid"][0] = XH2Of_initial ## mole fraction in the vapor
    df_results["XS_fluid"][0] = 0 ## vapor starts without sulfur
    df_results["phi_H2O"][0] = phi.phiH2O ## initial fugacity coefficient calculated at P_initial
    df_results["phi_H2S"][0] = phi.phiH2S ## initial fugacity coefficient calculated at P_initial
    df_results["phi_SO2"][0] = phi.phiSO2 ## initial fugacity coefficient calculated at P_initial
    df_results["S6+/ST"][0] = rs_melt_initial
    df_results["water_fugacity"][0] = fH2O_initial ## initial water fugacity in bar
    df_results["melt_fraction"][0] = 1 ## initial mass fraction of melt
    df_results["vapor_fraction"][0] = 0.0000 ## initial mass fraction of vapor
    df_results["crystal_fraction"][0] = 0.0000 ## initial mass fraction of crystal
    df_results["electron_balance"][0] = e_balance_initial ## redox budget relative to Fe3+ and S6+
    df_results["ferric"][0] = ferric_ratio_0 * melt_comp_initial.composition["FeOT"] / (55.845 + 15.999) ## initial moles of ferric iron in 100 g
    df_results["ferrous"][0] = (1 - ferric_ratio_0) * melt_comp_initial.composition["FeOT"] / (55.845 + 15.999) ## ferrous iron in 100 g
    df_results["ferric_ratio"][0] = ferric_ratio_0 ## Fe3+/FeT
    df_results["FeOT"][0] = melt_comp_initial.composition["FeOT"] ## in wt.%
    df_results["ferric_cr"][0] = 0 ## amount of Fe3+ taken by crystallization in moles in 100 g
    df_results["ferrous_cr"][0] = 0 ## amount of Fe2+ taken by crystallization in moles in 100 g
    df_results["FMQ"][0] = fo2_0.fmq() ## log10fO2 along FMQ buffer

    # degassing start
    m =l-0
    ## Decompression degassing
    for i in range(1, m):
        ## redefine the fugacity coefficient object at new pressure i, only related to P and T
        phi_volatiles = Fugacity(df_results["pressure"][i], temperature)
        df_results["phi_H2O"][i] = phi_volatiles.phiH2O
        df_results["phi_H2S"][i] = phi_volatiles.phiH2S
        df_results["phi_SO2"][i] = phi_volatiles.phiSO2

        # define melt composition, COH-only degassing, fo2, and sulfur partition coefficient objects with current pressure i, and melt fraction from previous step
        silicate_melt = MeltComposition(df_results["melt_fraction"][i - 1], choice)
        if COH_model == 1:
            coh_degas = VolatileCalc(tk, silicate_melt.composition["SiO2"])

        else:
            coh_degas = IaconoMarziano(df_results["pressure"][i], tk, silicate_melt.composition)

        fo2_degassing = OxygenFugacity(df_results["pressure"][i], tk, silicate_melt.composition)
        re = PartitionCoefficient(df_results["pressure"][i], tk, silicate_melt.composition, df_results["wH2O_melt"][i - 1],
                                  phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2)
        Fe2_cr = 0
        Fe3_cr = 0

        if fo2_tracker == 0: ## if fO2 tracker is disabled
            df_results["fO2"][i] = fo2_degassing.fo2(ferric_ratio_0) ## absolute log10fO2 is recalculated at each pressure using the same ferric/FeT ratio
            ## calculate three kds, SO2/ST in the vapor using the water fugacity from previous step, and fO2 from this pressure step.
            df_results["kd_RxnI"][i] = re.kd_rxn1(xh2o=df_results["XH2O_fluid"][i-1])
            df_results["kd_RxnII"][i] = re.kd_rxn2(fo2=10**(df_results["fO2"][i]))
            df_results["kd_RxnIa"][i] = re.kd_rxn1a(fo2=10**df_results["fO2"][i])
            df_results["SO2/ST"][i] = re.gas_quilibrium(fo2=10**df_results["fO2"][i], fh2o=df_results["water_fugacity"][i-1], phiso2=phi_volatiles.phiSO2,
                                                 phih2s=phi_volatiles.phiH2S)

            print(df_results["kd_RxnI"][i], df_results["kd_RxnII"][i], df_results["kd_RxnIa"][i])
            ## combined molar Kd weighed by SO2/ST in the vapor and S6+/ST in the melt
            df_results["kd_combined_molar"][i] = (df_results["SO2/ST"][i] * df_results["kd_RxnIa"][i] +
                                                 (1 - df_results["SO2/ST"][i]) * df_results["kd_RxnI"][i]) * (1 - rs_melt_initial) \
                                + rs_melt_initial * df_results["kd_RxnII"][i]

            ## S mole fraction in the vapor using combined molar Kd and mole fraction of S in the melt from previous step
            df_results["XS_fluid"][i] = df_results["XS_melt"][i-1] * df_results["kd_combined_molar"][i]
            df_results["XSO2_fluid"][i] = df_results["XS_fluid"][i] * df_results["SO2/ST"][i]
            df_results["XH2S_fluid"][i] = df_results["XS_fluid"][i] * (1 - df_results["SO2/ST"][i])
            print(df_results["XS_fluid"][i], df_results["kd_combined_molar"][i])

            if choice == 1: ## if crystallization is enabled
                ##With known S contents in the melt, using COH model and mass balance to solve for the CO2, H2O in the melt and in the vapor, and mass fractions of vapor, melt and crystal
                ## initial guess comes from the previous degassing step
                initial_guess = np.array(
                    [df_results["melt_fraction"][i-1], df_results["vapor_fraction"][i-1],
                     df_results["XH2O_fluid"][i-1], df_results["XCO2_fluid"][i-1], df_results["wH2O_melt"][i - 1],
                     df_results["wCO2_melt"][i - 1], df_results["crystal_fraction"][i-1]])
                root = coh_degas.coh_solubility(Pm=df_results["pressure"][i], h2o_guess=df_results["wH2O_melt"][i - 1]+0.3, co2_0= CO2_initial,
                                                h2o_0=H2O_initial, XS_fluid=df_results["XS_fluid"][i], rS_fluid=df_results["SO2/ST"][i],
                                                u0=initial_guess, choice=choice)



                if root.x[6] <= 0:
                    df_results["crystal_fraction"][i] = 0
                else:
                    df_results["crystal_fraction"][i] = root.x[6]

            else: ## if crystallization is disabled
                initial_guess = np.array(
                    [df_results["melt_fraction"][i-1], df_results["vapor_fraction"][i-1],
                    df_results["XH2O_fluid"][i-1], df_results["XCO2_fluid"][i-1], df_results["wH2O_melt"][i - 1],
                     df_results["wCO2_melt"][i - 1]])
                root = coh_degas.coh_solubility(h2o_guess=df_results["wH2O_melt"][i - 1],co2_0= CO2_initial,
                                                h2o_0=H2O_initial, XS_fluid=df_results["XS_fluid"][i], rS_fluid=df_results["SO2/ST"][i],
                                                u0=initial_guess, choice=choice)
                df_results["crystal_fraction"][i] = 0

            fm = root.x[0]
            fv = root.x[1]
            XH2O_f = root.x[2]
            XCO2_f = root.x[3]
            wtH2O_m = root.x[4]
            wtCO2_m = root.x[5]

            print(fv, fm, XH2O_f, XCO2_f, wtH2O_m, wtCO2_m)

            ## update melt composition after calculating the melt fraction fm
            melt_comp_updated = MeltComposition(melt_fraction=fm, choice=choice)
            ## update mass fractions of melt, vapor, CO2, H2O contents in the vapor and melt
            df_results["melt_fraction"][i] = fm
            df_results["vapor_fraction"][i] = fv
            df_results["XH2O_fluid"][i] = XH2O_f
            df_results["XCO2_fluid"][i] = XCO2_f
            df_results["wH2O_melt"][i] = wtH2O_m
            # if wtH2O_m <= df_results["wH2O_melt"][i-1]:
            #     df_results["wH2O_melt"][i] = wtH2O_m
            # else:
            #     df_results["wH2O_melt"][i] = df_results["wH2O_melt"][i-1]
            df_results["wCO2_melt"][i] = wtCO2_m

            fo2_degassing = OxygenFugacity(df_results["pressure"][i], tk, melt_comp_updated.composition)
            re_update = PartitionCoefficient(df_results["pressure"][i], tk, melt_comp_updated.composition, wtH2O_m,
                                             phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2)
            df_results["wS_fluid"][i] = 100 * (df_results["XSO2_fluid"][i] + df_results["XH2S_fluid"][i]) * 32.065 / (
                    XH2O_f * 18.015 + XCO2_f * 44.01 + df_results["XSO2_fluid"][i] * 64 + df_results["XH2S_fluid"][i] * 34)
            df_results["kd_combined_wt"][i] = df_results["wS_fluid"][i] * 10000 / df_results["wS_melt"][i - 1]
            # df_results["kd_combined_wt"][i] = 80

            if df_results["vapor_fraction"][i] + df_results["crystal_fraction"][i] > 0:
                df_results["DS_bulk"][i] = df_results["kd_combined_wt"][i] * df_results["vapor_fraction"][i] / (df_results["vapor_fraction"][i] + df_results["crystal_fraction"][i])
            else:
                df_results["DS_bulk"][i] = df_results["kd_combined_wt"][i]

            df_results["wS_melt"][i] = S_initial / (fm * (1 - df_results["DS_bulk"][i]) + df_results["DS_bulk"][i])
            # df_results["wS_melt"][i] = (S_initial - fv * df_results["wS_fluid"][i] *10000)/fm
            df_results["XS_melt"][i] = (df_results["wS_melt"][i] / (10000 * 32.065)) / (
                    re_update.ntot + df_results["wS_melt"][i] / (10000 * 32.065) + df_results["wH2O_melt"][i] / 18.015
                    + df_results["wCO2_melt"][i] / (10000 * 44.01))
            df_results["water_fugacity"][i] = df_results["XH2O_fluid"][i]*df_results["pressure"][i]*phi_volatiles.phiH2O
            # print(df_results["wS_melt"][i])
            if choice == 1:  ## redox budget of the crystals assuming crystals always take the same Fe3+/FeT as the melt in the previous step
                Fe3_cr = df_results["ferric_ratio"][i - 1] * (
                            df_results["FeOT"][i - 1] * df_results["melt_fraction"][i - 1] - fm *
                            melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
                Fe2_cr = (1 - df_results["ferric_ratio"][i - 1]) * \
                         (df_results["FeOT"][i - 1] * df_results["melt_fraction"][i - 1] - fm *
                          melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
                # e_FeO_cr = (1 - df_results["ferric_ratio"][i - 1]) * \
                #            (df_results["FeOT"][i-1] * df_results["melt_fraction"][i-1] - fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
            else:  ## if crystallization is disabled, =0
                Fe2_cr = 0
                Fe3_cr = 0

        else:## if fO2 tracker is enabled, fO2 change due to S degassing and crystallization of Fe3+ and Fe2+
            fo2_tr[0] = 10 ** (fo2_degassing.fo2(df_results["ferric_ratio"][i - 1]))  ## initial fO2 of current step using Fe3+/FeT ratio from the previous P step
            fH2O_initial = df_results["XH2O_fluid"][i - 1] * df_results["pressure"][i] * phi_volatiles.phiH2O  ## water fugacity using the XH2O from previous P step
            rs_fluid_initial = re.gas_quilibrium(fo2=fo2_tr[0], fh2o=fH2O_initial, phiso2=phi_volatiles.phiSO2,
                                                 phih2s=phi_volatiles.phiH2S)
            rs_f_tr[0] = rs_fluid_initial  ## initial XSO2/XST in the fluid using fO2 and fH2O from previous step
            XH2O_fluid_tr[0] = df_results["XH2O_fluid"][i - 1]
            XCO2_fluid_tr[0] = df_results["XCO2_fluid"][i - 1]
            XSO2_f_tr[0] = df_results["XSO2_fluid"][i - 1]
            XH2S_f_tr[0] = df_results["XH2S_fluid"][i - 1]
            wH2O_m_tr[0] = df_results["wH2O_melt"][i - 1]
            wCO2_m_tr[0] = df_results["wCO2_melt"][i - 1]
            wS_m_tr[0] = df_results["wS_melt"][i - 1]
            kd1_tr[0] = re.kd_rxn1(xh2o=XH2O_fluid_tr[0])
            kd2_tr[0] = re.kd_rxn2(fo2=fo2_tr[0])
            kd1a_tr[0] = re.kd_rxn1a(fo2=fo2_tr[0])
            ferric_ratio_tr[0] = df_results["ferric_ratio"][i - 1]
            rs_melt_initial = 10 ** (8 * np.log10(ferric_ratio_tr[0] / (1 - ferric_ratio_tr[0])) - 2863 / tk + 6.8)
            # rs_melt_initial = 10**(8 * np.log10(ferric_ratio_tr[0] / (1 - ferric_ratio_tr[0])) + 8.7436 * 1000000/(tk ** 2) - 27703/tk + 20.273)
            rs_melt_initial = rs_melt_initial / (1 + rs_melt_initial)
            rs_m_tr[0] = rs_melt_initial
            XS_m_tr[0] = (df_results["wS_melt"][i - 1] / (10000 * 32.065)) / (
                    re.ntot + df_results["wS_melt"][i - 1] / (10000 * 32.065) + re.nh + df_results["wCO2_melt"][i - 1] / (
                    10000 * 44.01))
            XS6_m_tr[0] = rs_melt_initial * XS_m_tr[0]
            XS2_m_tr[0] = (1 - rs_melt_initial) * XS_m_tr[0]

            if df_results["pressure"][i] >= 5:
                kd_combined_tr[0] = (rs_fluid_initial * kd1a_tr[0] + (1 - rs_fluid_initial) * kd1_tr[0]) * (1 - rs_melt_initial) \
                                + rs_melt_initial * kd2_tr[0]
            else:
                kd_combined_tr[0] = df_results["kd_combined_molar"][i-1]+20

            XS_f_tr[0] = XS_m_tr[0] * kd_combined_tr[0]
            melt_fraction_tr[0] = df_results["melt_fraction"][i - 1]
            vapor_fraction_tr[0] = df_results["vapor_fraction"][i - 1]
            crystal_fraction_tr[0] = df_results["crystal_fraction"][i - 1]
            a = [0]
            for j in range(1, def_variables.n):
                ## calculate initial solution for the system from the previous iteration

                XH2S_f_initial = kd_combined_tr[j - 1] * (1 - rs_f_tr[j - 1]) * XS_m_tr[j - 1]
                XSO2_f_initial = kd_combined_tr[j - 1] * rs_f_tr[j - 1] * XS_m_tr[j - 1]
                XH2O_f_initial = XH2O_fluid_tr[j - 1]
                XCO2_f_initial = 1 - XSO2_f_initial - XH2S_f_initial - XH2O_f_initial

                # recalculate H2O, CO2 in the melt, melt fraction, vapor fraction and crystal fraction after S degassing
                if choice == 1:
                    initial_guess = np.array(
                        [melt_fraction_tr[j - 1], vapor_fraction_tr[j - 1], XH2O_f_initial,
                         XCO2_f_initial, wH2O_m_tr[j-1], wCO2_m_tr[j-1], crystal_fraction_tr[j - 1]])
                    root = coh_degas.coh_solubility(Pm=df_results["pressure"][i], h2o_guess=df_results["wH2O_melt"][i - 1]+0.3, co2_0= CO2_initial,
                                                h2o_0=H2O_initial, XS_fluid=XSO2_f_initial + XH2S_f_initial, rS_fluid=rs_f_tr[j - 1],
                                                u0=initial_guess, choice=choice)

                    if root.x[6] <= 0:
                        crystal_fraction_tr[j] = 0
                    else:
                        crystal_fraction_tr[j] = root.x[6]

                else:
                    initial_guess = np.array(
                        [melt_fraction_tr[j - 1], vapor_fraction_tr[j - 1], XH2O_f_initial,
                         XCO2_f_initial, wH2O_m_tr[j-1], wCO2_m_tr[j-1]])
                    root = coh_degas.coh_solubility(Pm=df_results["pressure"][i], h2o_guess=df_results["wH2O_melt"][i - 1], co2_0=CO2_initial,
                                                h2o_0=H2O_initial, XS_fluid=XSO2_f_initial + XH2S_f_initial, rS_fluid=rs_f_tr[j - 1],
                                                u0=initial_guess, choice=choice)
                    crystal_fraction_tr[j] = 0
                fm = root.x[0]
                fv = root.x[1]
                XH2O_f = root.x[2]
                XCO2_f = root.x[3]
                wtH2O_m = root.x[4]
                wtCO2_m = root.x[5]
                ## update melt composition, fm, fv (and fc), H2O, CO2 in the melt with new values
                melt_comp_updated = MeltComposition(melt_fraction=fm, choice=choice)
                melt_fraction_tr[j] = fm
                vapor_fraction_tr[j] = fv
                XH2O_fluid_tr[j] = XH2O_f
                XCO2_fluid_tr[j] = XCO2_f
                wH2O_m_tr[j] = wtH2O_m
                # if wtH2O_m <= df_results["wH2O_melt"][i - 1]:
                #     df_results["wH2O_melt"][i] = wtH2O_m
                # else:
                #     df_results["wH2O_melt"][i] = df_results["wH2O_melt"][i - 1]-0.01
                wCO2_m_tr[j] = wtCO2_m
                fo2_degassing = OxygenFugacity(df_results["pressure"][i], tk, melt_comp_updated.composition)
                re_update = PartitionCoefficient(df_results["pressure"][i], tk, melt_comp_updated.composition, wtH2O_m,
                                                 phi_volatiles.phiH2O, phi_volatiles.phiH2S, phi_volatiles.phiSO2)
                if choice == 1: ## redox budget of the crystals assuming crystals always take the same Fe3+/FeT as the melt in the previous step
                    Fe3_cr = df_results["ferric_ratio"][i - 1] * (df_results["FeOT"][i-1] * df_results["melt_fraction"][i-1] - fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
                    Fe2_cr = (1 - df_results["ferric_ratio"][i - 1]) * \
                               (df_results["FeOT"][i-1] * df_results["melt_fraction"][i-1] - fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
                    # e_FeO_cr = (1 - df_results["ferric_ratio"][i - 1]) * \
                    #            (df_results["FeOT"][i-1] * df_results["melt_fraction"][i-1] - fm * melt_comp_updated.composition["FeOT"]) / (55.845 + 15.999)
                else:## if crystallization is disabled, =0
                    Fe2_cr = 0
                    Fe3_cr = 0
                    # e_FeO_cr = 0

                wS_f_tr[j] = 100 * (XSO2_f_initial + XH2S_f_initial) * 32.065 / (
                        XH2O_f * 18.015 + XCO2_f * 44.01 + XSO2_f_initial * 64 + XH2S_f_initial * 34)
                kdS_wt = wS_f_tr[j] * 10000 / df_results["wS_melt"][i - 1]

                if vapor_fraction_tr[j] + crystal_fraction_tr[j]>0:
                    DS_wt = kdS_wt * vapor_fraction_tr[j] / (vapor_fraction_tr[j] + crystal_fraction_tr[j])
                else:
                    DS_wt = kdS_wt
                wtS_m = S_initial / (fm * (1 - DS_wt) + DS_wt)
                # wtS_m = (S_initial - fv* wS_f_tr[j]*10000)/fm
                XS_m_tr[j] = (wtS_m / (10000 * 32.065)) / (
                        re_update.ntot + wtS_m / (10000 * 32.065) + 2* wtH2O_m / 18.015 + wtCO2_m / (10000 * 44.01))
                XS6_m_tr[j] = XS_m_tr[j] * rs_m_tr[j - 1]
                XS2_m_tr[j] = XS_m_tr[j] * (1 - rs_m_tr[j - 1])

                cohs = S_fO2(df_results["pressure"][i], tk, melt_comp_updated.composition)
                ## calculate the fO2 at the current iteration using Nash (Fe3+/FeT, S6+/ST), gas equilibrium (fO2, SO2/ST),
                # Kress and Carmicheal (Fe3+/FeT, fO2) and redox budget conservation
                ## Initial guess of S6+/ST, SO2/ST, Fe3+/FeT and fO2 are from previous iteration.
                low_l = 0.0001
                up_l = 0.9999

                if low_l <rs_m_tr[j - 1] < up_l:
                    if low_l < rs_f_tr[j - 1] < up_l:

                        initial_guess_2 = np.array(
                            [rs_m_tr[j - 1], rs_f_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                        [rs_m_n, rs_f_n, ferric_ratio_n, fo2_n] = \
                            cohs.cohs_solubility(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial, XS_m=XS_m_tr[j],
                                                 wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                                 phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                                 fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                                 u0=initial_guess_2)

                    elif rs_f_tr[j-1] > up_l:
                        initial_guess_2 = np.array([rs_m_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                        [rs_m_n, ferric_ratio_n, fo2_n] = \
                            cohs.cohs_so2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                             XS_m=XS_m_tr[j],
                                             wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                             phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                             fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                             u0=initial_guess_2)
                        rs_f_n = 1
                    else:
                        initial_guess_2 = np.array(
                            [rs_m_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]),
                             fo2_tr[j - 1]])
                        [rs_m_n, ferric_ratio_n, fo2_n] = \
                            cohs.cohs_h2s (fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                                 XS_m=XS_m_tr[j],
                                                 wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                                 phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                                 fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                                 u0=initial_guess_2)
                        rs_f_n = 0

                elif rs_m_tr[j - 1] > up_l:
                    if rs_f_tr[j-1] > up_l:
                        initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                        [ferric_ratio_n, fo2_n] = \
                            cohs.cohs_S6_so2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                             XS_m=XS_m_tr[j],
                                             wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                             phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                             fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                             u0=initial_guess_2)
                        rs_f_n = 1
                        rs_m_n = 1
                    elif low_l < rs_f_tr[j-1] < up_l:
                        initial_guess_2 = np.array(
                            [rs_f_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])
                        [rs_f_n, ferric_ratio_n, fo2_n] = \
                            cohs.cohs_solubility_S6(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                                    XS_m=XS_m_tr[j],
                                                    wS_f=wS_f_tr[j], wS_m=wS_m_tr[j],
                                                    wFeO_m=melt_comp_updated.composition["FeOT"],
                                                    phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                                    fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                                    u0=initial_guess_2)
                        rs_m_n = 1
                    else:
                        initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                        [ferric_ratio_n, fo2_n] = \
                            cohs.cohs_S6_so2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                             XS_m=XS_m_tr[j],
                                             wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                             phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                             fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                             u0=initial_guess_2)
                        rs_f_n = 0
                        rs_m_n = 1
                else:
                    if rs_f_tr[j-1] < low_l:
                        initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                        [ferric_ratio_n, fo2_n] = \
                            cohs.cohs_S2_h2s(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                             XS_m=XS_m_tr[j],
                                             wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                             phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                             fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                             u0=initial_guess_2)
                        rs_f_n = 0
                        rs_m_n = 0
                    elif low_l <rs_f_tr[j-1] < up_l:
                        initial_guess_2 = np.array(
                            [rs_f_tr[j - 1], ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])
                        [rs_f_n, ferric_ratio_n, fo2_n] = \
                            cohs.cohs_solubility_S2(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                                 XS_m=XS_m_tr[j],
                                                 wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                                 phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                                 fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                                 u0=initial_guess_2)
                        rs_m_n = 0
                    else:
                        initial_guess_2 = np.array([ferric_ratio_tr[j - 1] / (1 - ferric_ratio_tr[j - 1]), fo2_tr[j - 1]])

                        [ferric_ratio_n, fo2_n] = \
                            cohs.cohs_S2_h2s(fm=fm, fv=fv, XH2O_f=XH2O_f, XS_f=XSO2_f_initial + XH2S_f_initial,
                                             XS_m=XS_m_tr[j],
                                             wS_f=wS_f_tr[j], wS_m=wS_m_tr[j], wFeO_m=melt_comp_updated.composition["FeOT"],
                                             phi_so2=phi.phiSO2, phi_h2o=phi.phiH2O, phi_h2s=phi.phiH2S,
                                             fo2_cons=fo2_degassing.con, feo_cr_acc= df_results["ferrous_cr"][i-1], e_feo_cr=Fe2_cr, ebalance=e_balance_initial,
                                             u0=initial_guess_2)
                        rs_f_n = 1
                        rs_m_n = 0
                print(df_results["pressure"][i], rs_m_n, rs_f_n, ferric_ratio_n, np.log10(fo2_n), np.log10(fo2_tr[j - 1]),
                          wtH2O_m, fv, fm)

                ferric_ratio_tr[j] = ferric_ratio_n / (ferric_ratio_n + 1)
                fo2_tr[j] = fo2_n
                rs_f_tr[j] = rs_f_n
                rs_m_tr[j] = rs_m_n

                kd1_tr[j] = re_update.kd_rxn1(xh2o=XH2O_fluid_tr[j])
                # kd1_tr[j] = 50
                kd2_tr[j] = re_update.kd_rxn2(fo2=fo2_tr[j])
                kd1a_tr[j] = re_update.kd_rxn1a(fo2=fo2_tr[j])
                if df_results["pressure"][i] >= 5:
                    kd_combined_n = (kd1_tr[j] * (1 - rs_f_n) + kd1a_tr[j] * rs_f_n) * (1 - rs_m_n) + kd2_tr[j] * rs_m_n
                else:
                    kd_combined_n = df_results["kd_combined_molar"][i-1] + 20
                kd_combined_tr[j] = kd_combined_n
                XS_f_tr[j] = kd_combined_tr[j] * df_results["XS_melt"][i - 1]
                XSO2_f_tr[j] = XS_f_tr[j] * rs_f_tr[j]
                XH2S_f_tr[j] = XS_f_tr[j] * (1 - rs_f_tr[j])

                wS_f_tr[j] = 100 * (XSO2_f_tr[j] + XH2S_f_tr[j]) * 32.065 / (
                        XH2O_fluid_tr[j] * 18.015 + XCO2_fluid_tr[j] * 44.01 + XSO2_f_tr[j] * 64 + XH2S_f_tr[j] * 34)

                kdS_wt = wS_f_tr[j] * 10000 / df_results["wS_melt"][i - 1]

                if vapor_fraction_tr[j] + crystal_fraction_tr[j] > 0:
                    DS_wt = kdS_wt * vapor_fraction_tr[j] / (vapor_fraction_tr[j] + crystal_fraction_tr[j])
                else:
                    DS_wt = kdS_wt
                wtS_m = S_initial / (melt_fraction_tr[j] * (1 - DS_wt) + DS_wt)
                # wtS_m = (S_initial - fv * wS_f_tr[j] * 10000) / fm
                wS_m_tr[j] = wtS_m
                XS_m_tr[j] = (wtS_m / (10000 * 32.065)) / (
                        re_update.ntot + wtS_m / (10000 * 32.065) + 2 * wtH2O_m / 18.015 + wtCO2_m / (10000 * 44.01))
                XS6_m_tr[j] = XS_m_tr[j] * rs_m_tr[j]
                XS2_m_tr[j] = XS_m_tr[j] * (1 - rs_m_tr[j])

                # if abs(XS_m_tr[j]-XS_m_tr[j-1]) < 0.0001:
                if abs(np.log10(fo2_tr[j]) - np.log10(fo2_tr[j - 1])) < 0.01:
                    df_results["kd_RxnI"][i] = kd1_tr[j]
                    df_results["kd_RxnIa"][i] = kd1a_tr[j]
                    df_results["kd_RxnII"][i] = kd2_tr[j]
                    df_results["SO2/ST"][i] = rs_f_tr[j]
                    df_results["kd_combined_molar"][i] = kd_combined_tr[j]

                    df_results["kd_combined_wt"][i] = kdS_wt
                    df_results["DS_bulk"][i] = DS_wt
                    df_results["XS_melt"][i] = XS_m_tr[j]
                    df_results["wS_melt"][i] = wS_m_tr[j]
                    df_results["wS_fluid"][i] = wS_f_tr[j]
                    df_results["XS_fluid"][i] = XH2S_f_tr[j] + XSO2_f_tr[j]
                    df_results["XH2S_fluid"][i] = XH2S_f_tr[j]
                    df_results["XSO2_fluid"][i] = XSO2_f_tr[j]
                    df_results["SO2_fugacity"][i] = XSO2_f_tr[j] * phi_volatiles.phiSO2 * df_results["pressure"][i] * 10
                    df_results["H2S_fugacity"][i] = XH2S_f_tr[j] * phi_volatiles.phiH2S * df_results["pressure"][i] * 10
                    df_results["XH2O_fluid"][i] = XH2O_fluid_tr[j]
                    df_results["XCO2_fluid"][i] = XCO2_fluid_tr[j]
                    if wH2O_m_tr[j] <= df_results["wH2O_melt"][i - 1]:
                        df_results["wH2O_melt"][i] = wH2O_m_tr[j]
                    else:
                        df_results["wH2O_melt"][i] = df_results["wH2O_melt"][i - 1]
                    # df_results["wH2O_melt"][i] = wH2O_m_tr[j]
                    df_results["wCO2_melt"][i] = wCO2_m_tr[j]
                    df_results["melt_fraction"][i] = melt_fraction_tr[j]
                    df_results["vapor_fraction"][i] = vapor_fraction_tr[j]
                    df_results["crystal_fraction"][i] = crystal_fraction_tr[j]
                    df_results["water_fugacity"][i] = df_results["XH2O_fluid"][i] * df_results["pressure"][i] * 10 * phi_volatiles.phiH2O
                    df_results["fO2"][i] = np.log10(fo2_tr[j])
                    df_results["FMQ"][i] = fo2_degassing.fmq()
                    df_results["ferric_ratio"][i] = ferric_ratio_tr[j]
                    df_results["S6+/ST"][i] = rs_m_tr[j]
                    df_results["FeOT"][i] = melt_comp_updated.composition["FeOT"]
                    # df_results["Fe_cr"][i] = df_results["Fe_cr"][i-1]+e_FeO_cr

                    break
                elif j == def_variables.n:
                    df_results["kd_RxnI"][i] = kd1_tr[j]
                    df_results["kd_RxnIa"][i] = kd1a_tr[j]
                    df_results["kd_RxnII"][i] = kd2_tr[j]
                    df_results["SO2/ST"][i] = rs_f_tr[j]
                    df_results["kd_combined_molar"][i] = kd_combined_tr[j]
                    df_results["kd_combined_wt"][i] = kdS_wt
                    df_results["DS_bulk"][i] = DS_wt
                    df_results["XS_melt"][i] = XS_m_tr[j]
                    df_results["wS_melt"][i] = wS_m_tr[j]
                    df_results["wS_fluid"][i] = wS_f_tr[j]
                    df_results["XS_fluid"][i] = XH2S_f_tr[j] + XSO2_f_tr[j]
                    df_results["XH2S_fluid"][i] = XH2S_f_tr[j]
                    df_results["XSO2_fluid"][i] = XSO2_f_tr[j]
                    df_results["SO2_fugacity"][i] = XSO2_f_tr[j] * phi_volatiles.phiSO2 * df_results["pressure"][i] * 10
                    df_results["H2S_fugacity"][i] = XH2S_f_tr[j] * phi_volatiles.phiH2S * df_results["pressure"][i] * 10
                    df_results["XH2O_fluid"][i] = XH2O_fluid_tr[j]
                    df_results["XCO2_fluid"][i] = XCO2_fluid_tr[j]
                    df_results["wH2O_melt"][i] = wH2O_m_tr[j]
                    df_results["wCO2_melt"][i] = wCO2_m_tr[j]
                    df_results["melt_fraction"][i] = melt_fraction_tr[j]
                    df_results["vapor_fraction"][i] = vapor_fraction_tr[j]
                    df_results["crystal_fraction"][i] = crystal_fraction_tr[j]
                    df_results["water_fugacity"][i] = df_results["XH2O_fluid"][i] * df_results["pressure"][
                        i] * phi_volatiles.phiH2O
                    df_results["fO2"][i] = fo2_tr[j]
                    df_results["ferric_ratio"][i] = ferric_ratio_tr[j]
                    df_results["S6+/ST"][i] = rs_m_tr[j]
                    df_results["FeOT"][i] = melt_comp_updated.composition["FeOT"]
                    # df_results["Fe_cr"][i] = df_results["Fe_cr"][i-1]+e_FeO_cr
                    df_results["FMQ"][i] = fo2_degassing.fmq()
                    print("Calculation do not converge")

        melt_comp_updated = MeltComposition(df_results["melt_fraction"][i], choice)
        df_results["sulfate_m"][i] = df_results["melt_fraction"][i] * (df_results["wS_melt"][i] / 10000) * \
                                     df_results["S6+/ST"][i] / 32.065
        df_results["sulfide_m"][i] = df_results["melt_fraction"][i] * (df_results["wS_melt"][i] / 10000) * \
                                     (1 - df_results["S6+/ST"][i]) / 32.065
        df_results["SO2_f"][i] = df_results["vapor_fraction"][i] * df_results["wS_fluid"][i] * df_results["SO2/ST"][
            i] / 32.065
        df_results["H2S_f"][i] = df_results["vapor_fraction"][i] * df_results["wS_fluid"][i] * (
                1 - df_results["SO2/ST"][i]) / 32.065
        df_results["ferric"][i] = melt_comp_updated.composition["FeOT"] * df_results["melt_fraction"][i] * \
                                  df_results["ferric_ratio"][i] / (55.845 + 15.999)
        df_results["ferrous"][i] = melt_comp_updated.composition["FeOT"] * df_results["melt_fraction"][i] * \
                                   (1 - df_results["ferric_ratio"][i]) / (55.845 + 15.999)
        df_results["ferric_cr"][i] = df_results["ferric_cr"][i-1] + Fe3_cr
        df_results["ferrous_cr"][i] = df_results["ferrous_cr"][i-1] + Fe2_cr
        df_results["electron_balance"][i] = df_results["sulfide_m"][i] * 8 + df_results["H2S_f"][i] * 8 \
                                            + 2 * df_results["SO2_f"][i] + df_results["ferrous"][i] + df_results["ferrous_cr"][i]
    df_S_m[f"{m_run}"]=df_results["wS_melt"]
    df_CO2_m[f"{m_run}"] =df_results["wCO2_melt"]
df_S_m["average"]=df_S_m.mean(axis=1)
df_CO2_m["average"]=df_CO2_m.mean(axis=1)
df_S_m.to_csv ("Montecarlo_S_melt")
df_CO2_m.to_csv("Montecarlo_CO2_melt")

# df_results.to_csv("Fuego_results_ReMuth.csv")
# print(df_results["kd_RxnI"], df_results["kd_RxnIa"], df_results["kd_RxnII"], df_results["kd_combined_molar"],
#       df_results["kd_combined_wt"], df_results["DS_bulk"])
#

plt.figure(6)
plt.subplot(1, 2, 1)
plt.plot(df_S_m["average"][0:m], df_results["pressure"][0:m], "d")
plt.xlabel("S_melt (ppm)")
plt.ylabel("Pressure (MPa)")

plt.subplot(1, 2, 2)
plt.plot(df_S_m["average"][0:m], df_CO2_m["average"][0:m])
plt.plot(df["mi_S"], df["mi_CO2"], "o")
plt.xlabel("S_melt (ppm)")
plt.ylabel("CO2_melt (ppm)")
#
# plt.figure(7)
# plt.plot(df_results["fO2"][0:m]-df_results["FMQ"][0:m], df_results["pressure"][0:m], "v")
# plt.legend(["modeled dFMQ"])
# plt.xlabel("dFMQ")
# plt.ylabel("Pressure (MPa)")
#
# # plt.subplot(2, 2, 3)
# # plt.plot(df["water_melt"][1:m], df["carbon_dioxide_melt"][1:m], "o")
# # plt.plot(df_results["wH2O_melt"][1:m], df_results["wCO2_melt"][1:m], "v")
# # plt.xlabel("H2O_melt (wt.%)")
# # plt.ylabel("CO2_melt (ppm)")
# #
# # plt.subplot(2, 2, 4)
# m =m-1
# plt.figure(1)
# plt.plot(df_results["XH2O_fluid"][1:m], df_results["pressure"][1:m], "d")
# plt.plot(df_results["XCO2_fluid"][1:m], df_results["pressure"][1:m], "v")
# plt.plot(df_results["XH2S_fluid"][1:m], df_results["pressure"][1:m], "o")
# plt.plot(df_results["XSO2_fluid"][1:m], df_results["pressure"][1:m], "^")
# plt.legend(["XH2O_fluid", "XCO2_fluid", "XH2S_fluid", "XSO2_fluid"])
# plt.xlabel("mole fraction")
# plt.ylabel("Pressure (MPa)")
#
# plt.figure(2)
# plt.subplot(1, 2, 1)
# plt.plot(df_results["kd_RxnI"][0:m], df_results["pressure"][0:m], "d")
# plt.plot(df_results["kd_RxnIa"][0:m], df_results["pressure"][0:m], "v")
# plt.plot(df_results["kd_RxnII"][0:m], df_results["pressure"][0:m], "^")
# plt.plot(df_results["kd_combined_molar"][0:m], df_results["pressure"][0:m], "o")
# # plt.plot(df_results["kd_combined_wt"][0:m], df_results["pressure"][0:m], "o")
# plt.xlim([0, 150])
# plt.legend(["kdrxn1", " kdrxn1a", "kdrxn2", "kdcombined_molar"])
# plt.xlabel("partition coefficients")
# plt.ylabel("Pressure (MPa)")
#
# plt.subplot(1, 2, 2)
# plt.plot(df_results["melt_fraction"][0:m], df_results["pressure"][0:m], "d")
# plt.plot(df_results["vapor_fraction"][0:m], df_results["pressure"][0:m], "v")
# plt.plot(df_results["crystal_fraction"][0:m], df_results["pressure"][0:m], "o")
# plt.legend(["melt_fraction", "vapor_fraction", "crystal_fraction"])
# plt.xlabel("mass fraction")
# plt.ylabel("pressure (MPa)")
#
# plt.figure(4)
# plt.subplot(1, 2, 1)
# plt.plot(df_results["sulfate_m"][0:m], df_results["pressure"][0:m], "d")
# plt.plot(df_results["sulfide_m"][0:m], df_results["pressure"][0:m], "v")
# plt.legend(["sulfate_m", "sulfide_m"])
# plt.xlabel("moles (100g total)")
# plt.ylabel("Pressure (MPa)")
#
# plt.subplot(1, 2, 2)
# plt.plot(df_results["SO2_f"][0:m], df_results["pressure"][0:m], "^")
# plt.plot(df_results["H2S_f"][0:m], df_results["pressure"][0:m], "o")
# plt.legend(["SO2_f", "H2S_f"])
# plt.xlabel("moles (100g total)")
# plt.ylabel("Pressure (MPa)")
#
# plt.figure(5)
# # plt.subplot(1, 2, 1)
#
# plt.plot(df_results["ferric"][0:m], df_results["pressure"][0:m], "v")
# plt.plot(df_results["ferrous"][0:m], df_results["pressure"][0:m], "^")
# plt.plot(df_results["ferric_cr"][0:m], df_results["pressure"][0:m], "-")
# plt.plot(df_results["ferrous_cr"][0:m], df_results["pressure"][0:m], "+")
#
# plt.legend(["ferric", "ferrous", "ferric_crystallization", "ferrous_crystallization"])
# plt.xlabel("moles (100g total)")
# plt.ylabel("Pressure (MPa)")
#
# plt.figure(3)
# plt.subplot(2, 2, 1)
# plt.plot(df_results["ferric_ratio"][0:m], df_results["wH2O_melt"][0:m])
# # plt.plot(df_2["MK_Fe_r"], df_2["MK_H2O"], "d")
# # plt.plot(df_2["K_Fe_r"], df_2["K_H2O"], "s")
# # plt.legend(["Modeled", "Brounce et al 2017", "Moussallam et al., 2016"])
# plt.xlabel("Fe3+/FeT")
# plt.ylabel("H2O (wt.%)")
#
# plt.subplot(2, 2, 2)
# plt.plot(df_results["wS_melt"][0:m], df_results["wH2O_melt"][0:m])
# # plt.plot(df_2["MK_S"], df_2["MK_H2O"], "d")
# # plt.plot(df_2["K_S"], df_2["K_H2O"], "s")
# # plt.legend(["Modeled", "Brounce et al 2017", "Moussallam et al., 2016"])
# plt.xlabel("S (ppm)")
# plt.ylabel("H2O (wt.%)")
#
# plt.subplot(2, 2, 3)
# plt.plot(df_results["wS_melt"][0:m], df_results["S6+/ST"][0:m])
# # plt.plot(df_2["MK_S"], df_2["MK_S_r"], "d")
# # plt.legend(["Modeled", "Brounce et al 2017"])
# plt.ylabel("S6+/ST")
# plt.xlabel("S_melt (ppm)")
#
# plt.subplot(2, 2, 4)
# plt.plot(df_results["wS_melt"][0:m], df_results["ferric_ratio"][0:m])
# # plt.plot(df_2["MK_S"], df_2["MK_Fe_r"], "d")
# # plt.plot(df_2["K_S"], df_2["K_Fe_r"], "s")
# # plt.legend(["Modeled", "Brounce et al 2017", "Moussallam et al., 2016"])
# plt.xlabel("S_melt (ppm)")
# plt.ylabel("Fe3+/FeT")
#
# # plt.plot(df_results["water_fugacity"][0:m], df_results["pressure"][0:m], "d")
# # plt.plot(df_results["SO2_fugacity"][0:m], df_results["pressure"][0:m], "^")
# # plt.plot(df_results["H2S_fugacity"][0:m], df_results["pressure"][0:m], "o")
# # plt.xlabel("fugacity (bar)")
# # plt.ylabel("Pressure (MPa)")
# # plt.legend(["H2O", "SO2", "H2S"])
# plt.show()

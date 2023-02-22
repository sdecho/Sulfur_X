import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from oxygen_fugacity import OxygenFugacity
from fugacity import Fugacity
from sulfur_partition_coefficients import PartitionCoefficient
from Iacono_Marziano_COH import IaconoMarziano
from melt_composition import MeltComposition
from VC_COH import VolatileCalc
from newvariables import NewVariables
from degassingrun import COHS_degassing
from S_Fe import Sulfur_Iron
from SCSS_model import Sulfur_Saturation

# composition = {"SiO2": 50.42,
#                "Al2O3": 15.13,
#                "TiO2": 1.53,
#                "FeOT": 9.81,
#                "MgO": 7.76,
#                "CaO": 11.35,
#                "Na2O": 2.83,
#                "K2O": 0.14,
#                "P2O5": 0,
#                "MnO": 0,
#                }

# T = 1400
# P = 1000
# solubility = Sulfur_Saturation(P=P, T=T, sulfide_composition=sulfide, composition=composition, h2o=3.5, ferric_fe=0)
# scss = solubility.SCSS_smythe()
# scas = solubility.SCAS_Zajacz_Tsay()
# print(scss, scas)


# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

# _______________________________  USER INPUTS  ________________________________

# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
###############################################################################
# Basic input
###############################################################################
# Input of all the initial conditions
temperature = 1030  # temperature in C
# fO2 relative to FMQ buffer; if redox evolution is enabled, this is the initial fO2 at the initial P and T; if redox
# evolution is disabled, the degassing would be buffered at this delta_FMQ
delta_FMQ = 1.2
# initial H2O in wt.%
H2O_initial = 4.5
# initial CO2 in ppm
CO2_initial = 3300
# initial S in ppm
S_initial = 2650
# Crystallization or not? If 1, crystallization is enabled; if 0, crystallization disabled
choice = 1
# Which COH degassing model to use? 1 for VolatileCalc; 0 for Iacono_Marziano model
COH_model = 0
# Changing fO2 or not? If 1, fO2 changes by S degassing and S-Fe electron exchange; if 0, fO2 if buffered at the
# input delta_FMQ
fo2_tracker = 1
# Monte carlo simulation for S-degassing error estimate generated by errors in Kds can be performed by option.
# If 1, yes and please input the number of runs. This part can take time depending on the number of runs.
# If 0, the model stops after the first run.
monte_carlo = 0
if monte_carlo == 1:
    # If performing monte carlo simulation, please input the number of simulations.
    m_run = 500
else:
    m_run = 0
# name of the output file of monte_carlo simulations.
monte_name = "multiple_simulation_S.csv"
# The total steps of pressure from initial P to 1 bar. Please make sure the pressure step is small enough.
l = 1000
# total numbers of runs along degassing. This decides if the degassing model run till 1bar or not.
m = 1000
# sulfide composition (in wt.%), only relevant if SCSS is of interest;
# please change both here and in the degassingrun.py file
sulfide = {"Fe": 65.43,
            "Ni": 0,
            "Cu": 0,
            "O": 0,
            "S": 36.47
            }
###############################################################################
# Advanced parameters
# Please change the following parameters with caution and notify and justify
# your own choice if the results are used in publications.
###############################################################################
# Which S speciation model to use? 0 for Nash model; 1 for O'Neill and Mavrogenes (2022) model; 100 for Muth model;
# any other float number would be the modified version; the input is the last constant in the modified Muth model
S_Fe_choice = 0

# log10fO2 tolerance. The value of this number may cause individual outliners in the fO2 calcuation. However, if sigma
# is too small, it may cause the fault fO2 calculation
sigma = 0.005

# if crystallization is enabled, H2O-melt fraction relation is specified using H2O-K2O relation (K2O = a * H2O +b),
# assuming K2O is perfectly incompatible. The given a and b are based on H2O-K2O relation for Fuego magam from
# Lloyd et al. (2013). Both H2O and K2O are in wt.%. If crystallization is disabled, or running on magmas similar to Fuego,
# leave them unchanged.
slope_h2o = -0.713
constant_h2o = 3.689

###############################################################################
# Read-in CSV with MI data and model output
###############################################################################
# MI data and output csv file
# Read melt inclusion data that needs to be compared with the model results
# All the file have to be csv files with the format as in the example; file names need to end in ".csv"
mi_name = "Fuego.csv"
df = pd.read_csv(mi_name)
# name of the output csv file
output_name = f"Fuego_crystallization_{S_Fe_choice}.csv"

# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

# _______________________________  PERFORM DEGASSING________________________________

# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# The degassing run starts here.
# Define the inital variables
# Initiate composition class to calculate melt composition. If crystallization is enabled, melt composition is a
# function of degree of crystallization. If not, melt composition is constant
# Composition can be change in the file melt_composition
melt_comp_initial = MeltComposition(melt_fraction=1, choice=choice)

tk = float(temperature) + 273.15  # in kelvin
carbonate_initial = CO2_initial * 60.009 / 44.01

if COH_model == 1:  # if VolatileCalc is chosen
    VC_coh = VolatileCalc(TK=tk, sio2=melt_comp_initial.composition["SiO2"], a=slope_h2o, b=constant_h2o)
    [P_initial, WtH2O, PPMCO2, WtH2Om, WtOHm, XH2Of_initial, XCO2] = VC_coh.SatPress(WtH2O=H2O_initial,
                                                                                     PPMCO2=CO2_initial)
    print(
        f" The initial vapor saturation pressure is {P_initial} bar, and the initial vapor concentration is XH2O = {XH2Of_initial}"
        f" and XCO2 = {1 - XH2Of_initial}")
else:  # if IaconoMarziano model is chosen
    coh = IaconoMarziano(pressure=400, temperature_k=tk, composition=melt_comp_initial.composition,
                         a=slope_h2o, b=constant_h2o)
    [P_initial, XH2Of_initial] = coh.saturation_pressure(carbonate_initial, H2O_initial)  # P_initial in bar
    print(
        f" The initial vapor saturation pressure is {P_initial} bar, and the initial vapor concentration is XH2O = {XH2Of_initial}"
        f" and XCO2 = {1 - XH2Of_initial}.")
    coh = IaconoMarziano(pressure=P_initial / 10, temperature_k=tk, composition=melt_comp_initial.composition,
                         a=slope_h2o, b=constant_h2o)

# Define the variables
def_variables = NewVariables(P_initial, l)
my_data = def_variables.results_dic()
# Dataframe to store the output results
df_results = pd.DataFrame(data=my_data)
# This dataframe will only be used if the monte-carlo simulation option is enabled.
df_results_m = pd.DataFrame(data=my_data)

# initiate OxygenFugacity module to calculate fO2 with given ferric iron ratio or vice versa
fo2_0 = OxygenFugacity(P_initial / 10, tk, melt_comp_initial.composition)
# float(input("What is the initial Ferric/total Fe ratio in the melt?")) # Fe3+/FeT
ferric_ratio_0 = fo2_0.fe_ratio((fo2_0.fmq() + delta_FMQ))
print(ferric_ratio_0)
# initiate Fugacity class to calculate fugacity coefficients for H2O, SO2, H2S
phi = Fugacity(P_initial / 10, temperature)
# initiate PartitionCoefficient class to calculate sulfur partition coefficients
re = PartitionCoefficient(P_initial / 10, tk, melt_comp_initial.composition, H2O_initial, phi.phiH2O, phi.phiH2S,
                          phi.phiSO2, monte=0)
# initiate the solubility module to calculate the SCSS and SCAS
solubility = Sulfur_Saturation(P=P_initial/10, T=temperature, sulfide_composition=sulfide,
                               composition=melt_comp_initial.composition, h2o=H2O_initial, ferric_fe=ferric_ratio_0)
# Initial value of different variables
# Initial sulfur mole fraction in the melt
XS_initial = (S_initial / (10000 * 32.065)) / (
        re.ntot + S_initial / (10000 * 32.065) + re.nh + CO2_initial / (10000 * 44.01))
# initial water fugacity calculated by initial water contents in the vapor
fH2O_initial = XH2Of_initial * P_initial * phi.phiH2O
fH2_initial = re.hydrogen_equilibrium(fh2o=fH2O_initial,fo2=10**(fo2_0.fo2(ferric_ratio_0)))
# initiate Sulfur-Iron to calculate the sulfur speciation with Fe3+/FeT and T input
rs_melt = Sulfur_Iron(ferric_iron=ferric_ratio_0, temperature=temperature, model_choice=S_Fe_choice,
                      composition=melt_comp_initial.composition, o2=fo2_0.fmq() + delta_FMQ)
rs_melt_initial = rs_melt.sulfate

# moles of RB (redox budget) relative to S6+ and Fe3+ in 100 g of starting material
e_balance_initial = (S_initial / 10000) * (1 - rs_melt_initial) * 8 / 32.065 \
                    + (1 - ferric_ratio_0) * melt_comp_initial.composition["FeOT"] / (55.845 + 15.999)
print(rs_melt_initial)
# define all the initial values in the results DF.
df_results["wS_melt"][0] = S_initial  # in ppm
df_results["wH2O_melt"][0] = H2O_initial  # in wt.%
df_results["wCO2_melt"][0] = CO2_initial  # in ppm
df_results["XS_melt"][0] = XS_initial  # mole fraction in the melt
df_results["fO2"][0] = fo2_0.fo2(ferric_ratio_0)  # log10fO2
df_results["XCO2_fluid"][0] = 1 - XH2Of_initial  # mole fraction, in the vapor
df_results["XH2O_fluid"][0] = XH2Of_initial  # mole fraction in the vapor
df_results["XS_fluid"][0] = 0  # vapor starts without sulfur
df_results["phi_H2O"][0] = phi.phiH2O  # initial fugacity coefficient calculated at P_initial
df_results["phi_H2S"][0] = phi.phiH2S  # initial fugacity coefficient calculated at P_initial
df_results["phi_SO2"][0] = phi.phiSO2  # initial fugacity coefficient calculated at P_initial
df_results["S6+/ST"][0] = rs_melt_initial  # initial S6+/ST
df_results["water_fugacity"][0] = fH2O_initial  # initial water fugacity in bar
df_results["melt_fraction"][0] = 1  # initial mass fraction of melt
df_results["vapor_fraction"][0] = 0.0000  # initial mass fraction of vapor
df_results["crystal_fraction"][0] = 0.0000  # initial mass fraction of crystal
df_results["electron_balance"][0] = e_balance_initial  # redox budget relative to Fe3+ and S6+
# initial moles of ferric iron in 100 g
df_results["ferric"][0] = ferric_ratio_0 * melt_comp_initial.composition["FeOT"] / (55.845 + 15.999)
# initial moles of ferrous iron in 100 g
df_results["ferrous"][0] = (1 - ferric_ratio_0) * melt_comp_initial.composition["FeOT"] / (55.845 + 15.999)
df_results["ferric_ratio"][0] = ferric_ratio_0  # initial Fe3+/FeT
df_results["FeOT"][0] = melt_comp_initial.composition["FeOT"]  # in wt.%
df_results["ferric_cr"][0] = 0  # amount of Fe3+ taken by crystallization in moles in 100 g
df_results["ferrous_cr"][0] = 0  # amount of Fe2+ taken by crystallization in moles in 100 g
df_results["FMQ"][0] = fo2_0.fmq()  # initial log10fO2 along FMQ buffer
df_results["SCSS"][0] = solubility.SCSS_smythe()
df_results["SCAS"][0] = solubility.SCAS_Zajacz_Tsay()
df_results["fH2"][0] = fH2_initial

df_results_m = df_results

# degassing start
for i in range(1, m):
    degas = COHS_degassing(pressure=df_results["pressure"][i], temperature=temperature, COH_model=COH_model,
                           xlt_choice=choice, S_Fe_choice=S_Fe_choice, H2O_initial=H2O_initial, CO2_initial=CO2_initial,
                           S_initial=S_initial, a=slope_h2o, b=constant_h2o, monte_c=0)
    if fo2_tracker == 1:
        df_results.iloc[i] = degas.degassing_redox(df_results=df_results, index=i, e_balance_initial=e_balance_initial,
                                                   sigma=sigma)
    else:
        df_results.iloc[i] = degas.degassing_noredox(df_results=df_results, index=i, delta_FMQ=delta_FMQ)

# save outputs to a csv file
df_results.to_csv(output_name)

###############################################################################
# Plot the results
###############################################################################
# plot results

# This figure is similar to Figure 6 in the text
plt.figure(20)
plt.subplot(1,2,1)
plt.subplot(1, 2, 1)
plt.plot(df_results["wS_melt"][0:m], df_results["pressure"][0:m], linestyle="-", linewidth=5,
         color="blue")
plt.xlabel("S_melt (ppm)")
plt.ylabel("Pressure (MPa)")
plt.xlim([0,3000])
plt.ylim([0, 600])
plt.subplot(1,2,2)
plt.plot(df_results["fO2"][13:m] - df_results["FMQ"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=5,
         color="blue")
plt.legend(["modeled dFMQ"])
plt.xlabel("dFMQ")
plt.ylabel("Pressure (MPa)")
plt.xlim([0, 2])
plt.ylim([0, 600])



# This is similar to the Figure 7 in the main text
plt.figure(7)
plt.subplot(3, 2, 1)
plt.plot(df_results["kd_RxnI"][0:m], df_results["pressure"][0:m], linestyle="--", linewidth=3)
# plt.plot(df_results["kd_RxnIa"][13:m], df_results["pressure"][13:m], linestyle="-.", linewidth=3)
plt.plot(df_results["kd_RxnII"][13:m], df_results["pressure"][13:m], linestyle=":", linewidth=3)
plt.plot(df_results["kd_combined_molar"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=5, color="grey")
# plt.plot(df_results["kd_combined_wt"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=5, color="yellow")
plt.xlim([0, 300])
plt.ylim([0, 400])
plt.legend(["kdrxn1", "kdrxn2", "kdcombined_molar"])
plt.xlabel("partition coefficients")
plt.ylabel("Pressure (MPa)")

plt.subplot(3, 2, 2)
plt.plot(np.log10(df_results["kd_combined_wt"][1:m]), df_results["pressure"][1:m], linestyle="-", linewidth=5, color="yellow")
plt.legend(["This model"])
plt.xlabel("Kd_wt")
plt.ylabel("Pressure (MPa)")
plt.xlim([0, 4])
plt.ylim([0, 700])

plt.subplot(3, 2, 3)
plt.plot(df_results["ferric_ratio"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=3, color="grey")
plt.plot(df_results["SO2/ST"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=3, color="blue")
plt.plot(df_results["S6+/ST"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=3, color="red")
plt.xlim([0.2, 1])
plt.ylim([0, 700])
plt.legend(["Fe3+/FeT", "SO2/STV", "S6+/ST"])
plt.xlabel("ratios")
plt.ylabel("Pressure (MPa)")

plt.subplot(3, 2, 4)
plt.plot(df_results["fO2"][13:m] - df_results["FMQ"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=5,
         color="yellow")
plt.legend(["modeled dFMQ"])
plt.xlabel("dFMQ")
plt.ylabel("Pressure (MPa)")
plt.xlim([0, 2])
plt.ylim([0, 700])

plt.subplot(3, 2, 5)
plt.plot(np.log10(df_results["SO2_fugacity"][13:m]), df_results["pressure"][13:m], linestyle="-", linewidth=3,
         color="blue")
plt.plot(np.log10(df_results["H2S_fugacity"][13:m]), df_results["pressure"][13:m], linestyle="--", linewidth=2,
         color="blue")
plt.plot(np.log10(df_results["water_fugacity"][13:m]), df_results["pressure"][13:m], linestyle=":", linewidth=4,
         color="blue")
plt.legend(["fSO2", "fH2S", "fH2O"])
plt.xlabel("fugacity (bar)")
plt.ylabel("Pressure (MPa)")
plt.xlim([-2, 4])
plt.ylim([0, 700])

plt.subplot(3, 2, 6)
plt.plot(df_results["sulfate_m"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=3, color="red")
plt.plot(df_results["sulfide_m"][13:m], df_results["pressure"][13:m], linestyle="--", linewidth=2, color="red")
plt.plot(df_results["SO2_f"][13:m], df_results["pressure"][13:m], linestyle="-", linewidth=3, color="blue")
plt.plot(df_results["H2S_f"][13:m], df_results["pressure"][13:m], linestyle="--", linewidth=2, color="blue")
plt.legend(["S6+", "S2-", "SO2", "H2S"])
plt.xlim([0, 0.01])
plt.ylim([0, 700])
plt.xlabel("moles")
plt.ylabel("Pressure (MPa)")

# This figure is similar to Figure 8 in the text
plt.figure(8)
plt.subplot(1, 2, 1)
plt.plot(df_results["wS_melt"][0:m], df_results["pressure"][0:m], linestyle='-', linewidth=2)
plt.xlabel("S_melt (ppm)")
plt.ylabel("Pressure (MPa)")
# plt.xlim([200, 1600])
# plt.ylim([0, 700])

plt.subplot(1, 2, 2)
plt.plot(np.log10(df_results["kd_RxnI"][0:m]), df_results["pressure"][0:m], linestyle="--", linewidth=3)
# plt.plot(np.log10(df_results["kd_RxnIa"][0:m]), df_results["pressure"][0:m], linestyle="-.", linewidth=3)
plt.plot(np.log10(df_results["kd_RxnII"][0:m]), df_results["pressure"][0:m], linestyle=":", linewidth=3)
plt.plot(np.log10(df_results["kd_combined_molar"][0:m]), df_results["pressure"][0:m], linestyle="-", linewidth=5,
         color="grey")
plt.plot(np.log10(df_results["kd_combined_wt"][0:m]), df_results["pressure"][0:m], linestyle="-", linewidth=5,
         color="black")
# plt.xlim([-2, 5])
# plt.ylim([0, 700])
plt.legend(["kdrxn1", "kdrxn2", "kdcombined_molar", "kdcombined_wt"])
plt.xlabel("partition coefficients (log10)")
plt.ylabel("Pressure (MPa)")

# This figure is similar to Figure 8 in the text
# plt.figure(9)
# plt.subplot(2, 2, 1)
# plt.plot(df_results["wS_melt"][0:m], df_results["wH2O_melt"][0:m], linestyle='-', linewidth=2)
# plt.plot(df["mi_S"][0:30], df["mi_H2O"][0:30], markerfacecolor='green', marker='d', markeredgecolor="black",
#          markersize=10, linestyle="none")
# plt.plot(df["mi_S"][31:], df["mi_H2O"][31:], markerfacecolor='magenta', marker='s', markeredgecolor="black",
#          markersize=10, linestyle="none")
# plt.legend(["Modeled", "Brounce et al 2017", "Moussallam et al., 2016"])
# plt.xlabel("S (ppm)")
# plt.ylabel("H2O (wt.%)")
# plt.xlim([0, 1600])
# plt.ylim([0, 1])
#
# plt.subplot(2, 2, 2)
# plt.plot(df_results["wS_melt"][1:m], df_results["water_fugacity"][1:m] / 100, linestyle='-', linewidth=2)
# plt.plot(df_results["wS_melt"][1:m], df_results["SO2/ST"][1:m], linestyle='-', linewidth=2)
# plt.ylim([0, 2])
# plt.xlim([0, 1600])
# plt.legend(["fH2O", "SO2/ST"])
#
# plt.subplot(2, 2, 3)
# plt.plot(df_results["wS_melt"][0:m], df_results["S6+/ST"][0:m], linestyle='-', linewidth=2)
# plt.plot(df["mi_S"][0:30], df["mi_S_r"][0:30], markerfacecolor='green', marker='d', markeredgecolor="black",
#          markersize=10, linestyle="none")
# plt.legend(["Modeled", "Brounce et al 2017"])
# plt.ylabel("S6+/ST")
# plt.xlabel("S_melt (ppm)")
# plt.xlim([0, 1600])
# plt.ylim([-0.01, 0.1])
#
# plt.subplot(2, 2, 4)
# plt.plot(df_results["wS_melt"][0:m], df_results["ferric_ratio"][0:m], linestyle='-', linewidth=2)
# plt.plot(df["mi_S"][0:30], df["mi_Fe_r"][0:30], markerfacecolor='green', marker='d', markeredgecolor="black",
#          markersize=10, linestyle="none")
# plt.plot(df["mi_S"][30:], df["mi_Fe_r"][30:], markerfacecolor='magenta', marker='s', markeredgecolor="black",
#          markersize=10, linestyle="none")
# plt.legend(["Modeled", "Brounce et al 2017", "Moussallam et al., 2016"])
# plt.xlabel("S_melt (ppm)")
# plt.ylabel("Fe3+/FeT")
# plt.xlim([0, 1600])
# plt.ylim([0.12, 0.22])

# plt.show()

##############################################################################
# Monte-Carlo Simulation
##############################################################################
# df_S_m = df_results[["pressure"]].copy()
# print("Montecarlo simulation")
# if monte_carlo == 1:
#     for k in range(0, m_run):
#         for i in range(1, m):
#             degas = COHS_degassing(pressure=df_results_m["pressure"][i], temperature=temperature, COH_model=COH_model,
#                                    xlt_choice=choice, S_Fe_choice=S_Fe_choice, H2O_initial=H2O_initial,
#                                    CO2_initial=CO2_initial,
#                                    S_initial=S_initial, a=slope_h2o, b=constant_h2o, monte_c=1)
#             if fo2_tracker == 1:
#                 df_results_m.iloc[i] = degas.degassing_redox(df_results=df_results_m, index=i,
#                                                              e_balance_initial=e_balance_initial,
#                                                              sigma=sigma)
#             else:
#                 df_results_m.iloc[i] = degas.degassing_noredox(df_results=df_results_m, index=i, delta_FMQ=delta_FMQ)
#         df_S_m.insert(k + 1, k, df_results_m["wS_melt"])
# S_only = df_S_m.iloc[:, 1:m_run]
# df_S_m.insert(m_run + 1, "mean", S_only.mean(1))
# df_S_m.insert(m_run + 2, "std", S_only.std(1))
# df_S_m.insert(m_run + 3, "variance", S_only.var(1))
df_S_m = pd.read_csv("multiple_simulation_S_1_6000.csv")
# print(df_S_m.iloc[:, 1])
plt.figure(84)
plt.subplot(1, 2, 1)
plt.plot(df_results["wS_melt"][0:m], df_results["wCO2_melt"][0:m])
# for i in range(1, m_run):
#     plt.plot(df_S_m.iloc[:, i],df_results["wCO2_melt"])
# plt.plot(df["mi_S"], df["mi_CO2"], "o")
plt.legend(["Sulfur_X", " MI"])
plt.xlabel("S_melt (ppm)")
plt.ylabel("CO2_melt (ppm)")
plt.xlim([0,3000])
plt.ylim([0, 6000])
plt.subplot(1, 2, 2)
plt.plot(df_results["wS_melt"][0:m], df_results["wH2O_melt"][0:m])
# for i in range(1, m_run):
#     plt.plot(df_S_m.iloc[:, i],df_results["wH2O_melt"])
# plt.plot(df["mi_S"], df["mi_H2O"], "o")
# plt.legend(["Sulfur_X", " MI"])
plt.xlabel("S_melt (ppm)")
plt.ylabel("H2O_melt (wt.%)")
plt.xlim([0, 3000])
plt.ylim([0, 5])
plt.show()
# df_S_m.to_csv(monte_name)

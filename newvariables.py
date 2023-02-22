import numpy as np


class NewVariables:
    """
    P_initial: Initial pressure, in bar
    l: total pressure steps: int
    """
    def __init__(self, P_initial, l):
        self.n = 1000
        self.P = P_initial
        self.step = l

    def iteration_v(self, XH2Of_initial, ferric_ratio_0, H2O_initial, CO2_initial, S_initial):
        """
        XH2Of_initial: initial water fugacity in bar
        ferric_ratio_0: initial Fe3+/FeT (molar ratio)
        H2O_initial: initial H2O concentration in wt.%
        CO2_initial: initial CO2 concentration in wt.%
        S_initial: initial S concentration in ppm
        """
        fo2_tr = [0.5] * self.n
        XH2O_fluid_tr = [XH2Of_initial] * self.n
        XCO2_fluid_tr = [0.5] * self.n
        XSO2_f_tr = [0.5] * self.n
        XH2S_f_tr = [0.5] * self.n
        XS_f_tr = [0.5] * self.n
        ferric_ratio_tr = [ferric_ratio_0] * self.n
        wS_f_tr = [0.5] * self.n
        XS6_m_tr = [0.5] * self.n
        XS2_m_tr = [0.5] * self.n
        XS_m_tr = [0.1] * self.n
        wH2O_m_tr = [H2O_initial] * self.n
        wCO2_m_tr = [CO2_initial] * self.n
        wS_m_tr = [S_initial] * self.n
        kd1_tr = [0.5] * self.n
        kd2_tr = [0.5] * self.n
        kd1a_tr = [0.5] * self.n
        kd_combined_tr = [0.5] * self.n
        rs_m_tr = [0.5] * self.n
        rs_f_tr = [0.5] * self.n
        melt_fraction_tr = [0.8] * self.n
        crystal_fraction_tr = [0.1] * self.n
        vapor_fraction_tr = [0.1] * self.n
        return fo2_tr, XH2O_fluid_tr, XCO2_fluid_tr, XSO2_f_tr, XH2S_f_tr, XS_f_tr, ferric_ratio_tr, wS_f_tr, XS6_m_tr, \
               XS2_m_tr, XS_m_tr, wH2O_m_tr, wCO2_m_tr, wS_m_tr, kd1_tr, kd2_tr, kd1a_tr, kd_combined_tr, rs_m_tr, \
               rs_f_tr, melt_fraction_tr, crystal_fraction_tr, vapor_fraction_tr

    def results_dic(self):
        pressure = np.linspace(self.P / 10, 1, self.step)
        empty_list = np.linspace(0, 0, len(pressure))
        my_data = {"pressure": pressure,
                   "fO2": empty_list,
                   "wS_melt": empty_list,
                   "wH2O_melt": empty_list,
                   "wCO2_melt": empty_list,
                   "XS_melt": empty_list,
                   "phi_H2O": empty_list,
                   "XS_fluid": empty_list,
                   "XH2O_fluid": empty_list,
                   "XCO2_fluid": empty_list,
                   "XSO2_fluid": empty_list,
                   "XH2S_fluid": empty_list,
                   "phi_H2S": empty_list,
                   "phi_SO2": empty_list,
                   "wS_fluid": empty_list,
                   "melt_fraction": empty_list,
                   "vapor_fraction": empty_list,
                   "crystal_fraction": empty_list,
                   "DS_bulk": empty_list,
                   "kd_combined_wt": empty_list,
                   "kd_combined_molar": empty_list,
                   "kd_RxnI": empty_list,
                   "kd_RxnIa": empty_list,
                   "kd_RxnII": empty_list,
                   "SO2/ST": empty_list,
                   "S6+/ST": empty_list,
                   "water_fugacity": empty_list,
                   "SO2_fugacity": empty_list,
                   "H2S_fugacity": empty_list,
                   "electron_balance": empty_list,
                   "sulfate_m": empty_list,
                   "sulfide_m": empty_list,
                   "SO2_f": empty_list,
                   "H2S_f": empty_list,
                   "ferric": empty_list,
                   "ferrous": empty_list,
                   "ferric_ratio": empty_list,
                   "FeOT": empty_list,
                   "Fe_cr": empty_list,
                   "ferric_cr": empty_list,
                   "ferrous_cr": empty_list,
                   "FMQ": empty_list,
                   "SCSS":empty_list,
                   "SCAS":empty_list,
                   "fH2":empty_list,
                   }
        return my_data


class isotope:
    """
    T: in C
    f_s: the fraction of the S remaining in the melt at current degassing step
    f_sulfate: S6+/St in the melt
    f_so2: xso2/(xso2+xh2s) in the vapor phase
    d34s_ini: the initial δ34S in the melt
    """
    def __init__(self, Tc, f_sulfate, f_so2):
        self.t = Tc
        self.sulfate = f_sulfate
        self.so2 = f_so2 
        self.alpha_gas_melt = self.fractionation_factor()
        
        
    
    def fractionation_factor (self):
        #1000lna_so2_h2s and 1000lna_h2s_sulfide by Taylor 1986
        lna_so2_h2s = -0.42*(1000/self.t)**3+4.367*(1000/self.t)**2-0.105*1000/self.t-0.41
        lna_h2s_sulfide = 1.1*(1000/self.t)**2-0.19

        #1000lna_sulfate_sulfide and 1000lna_sulfate_h2s by Miyoshi et al.,1984
        lna_sulfate_sulfide = 7.4*(1000/self.t)**2-0.19
        lna_sulfate_h2s = 6.5*(1000/self.t)**2

        lna_so2_sulfate = lna_so2_h2s-lna_sulfate_h2s
        lna_so2_sulfide = lna_so2_sulfate+lna_sulfate_sulfide

        lna_gas_melt = self.so2*lna_so2_h2s-self.sulfate*lna_sulfate_sulfide+lna_h2s_sulfide
        return lna_gas_melt
    
 



#CPA e CPCA
import numpy as np

#Constantes
par={"R":8.3144598, "Nav":6.02214076e23, "vk":[0.48, 1.574, -0.176]} 

def alfaT(T, par):
    Tr = T / par["Tc"]
    return (1 + par["c1"] * (1 - np.sqrt(Tr)))**2

def Z(rho, T, x, par):
    return 1.51

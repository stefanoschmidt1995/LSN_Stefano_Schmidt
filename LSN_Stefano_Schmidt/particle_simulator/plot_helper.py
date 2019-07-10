#*************************************************************************
#STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
#EXERCISE 4/7 - HELPER CODE 
#*************************************************************************
import numpy as np
import matplotlib.pyplot as plt

def make_plots(gas_type = None, sigma = None, epsilon_kb = None, m =None):
    """Makes plots against time of various phyisical quantities: E_tot, E_kin, E_tot, P, T.
    It reads from file /out/avg_phasetype.dat with header
        [rho, [E_tot}, [E_kin], [E_pot], [virial], [Temp], [P], moves]
        [0,    1,2        3,4      5,6      7,8    9,10  11,12    13 ]
    Units: [sigma] = nm; [epsilon_kb] = K; [m]=amu
    """
    gas_label=""
    if gas_type is not None:
        gas_label = str(gas_type)
    phys_cond = ["solid","liquid", "gas"]
        #dealing with plots
    fig, energy_plots = plt.subplots(1,3, figsize=(15,6))
    fig2, PT_plots = plt.subplots(1,2,figsize=(15,6))
    fig.suptitle("Simulation of "+gas_label, size=20)
    
    k_B = 1.380#e-23
    for cond in phys_cond:
        data = np.loadtxt("./out/avg_"+cond+".dat", skiprows=1)
                #setting SI units (only values are considered; exponents are set by hands)
        if sigma is not None and epsilon_kb is not None:
            #1 amu = 1.66e-27 kg
            data[:,1] *= 1./sigma**3
            if m is not None: #if m is None time is thought to have the meaning of montecarlo step and thus is dimensionless
                data[:,13] *= sigma*np.sqrt(m*1.66/epsilon_kb*k_B)*(10) #time (ps)
            data[:,9:10] *= epsilon_kb #temperature (K)
            data[:,11:12] *= epsilon_kb*k_B/sigma**3 *10**-4 #pressure (Pa)
            data[:,1:2] *= (epsilon_kb*k_B)*6.022*10**-3 #E_tot (kJ/mol)
            data[:,3:4] *=(epsilon_kb*k_B)*6.022*10**-3 #E_kin (J/mol)
            data[:,5:6] *=(epsilon_kb*k_B)*6.022*10**-3 #E_pot (J/mol)
        energy_plots[0].errorbar(data[:,13], data[:,1], yerr= data[:,2],  fmt='-o',label="E_tot "+cond, markersize=2)
        energy_plots[1].errorbar(data[:,13], data[:,3], yerr= data[:,4],  fmt='-o',label="E_kin "+cond, markersize=2)
        energy_plots[2].errorbar(data[:,13], data[:,5], yerr= data[:,6],  fmt='-o',label="E_pot "+cond, markersize=2)
        PT_plots[0].errorbar(data[:,13], data[:,9], yerr= data[:,10],  fmt='-o',label="T "+cond, markersize=2)
        PT_plots[1].errorbar(data[:,13], data[:,11], yerr= data[:,12],  fmt='-o',label="P "+cond, markersize=2)

    for plot in energy_plots:
        if m is not None:
            plot.set_xlabel("Time ($10^{-12}$ s)")
        else:
            plot.set_xlabel("Montecarlo steps")
        plot.set_ylabel("Energy per particle (kJ/mol)")
        plot.legend()
    energy_plots[0].set_title("Total Energy "+gas_label)
    energy_plots[1].set_title("Kinetic Energy "+gas_label)
    energy_plots[2].set_title("Potential Energy "+gas_label)
    if m is not None:
        PT_plots[0].set_xlabel("Time ($10^{-12}$ s)")
        PT_plots[1].set_xlabel("Time ($10^{-12}$ s)")
    else:
        PT_plots[0].set_xlabel("Montecarlo steps")
        PT_plots[1].set_xlabel("Montecarlo steps")
    PT_plots[0].set_ylabel("Temperature (K)")
    PT_plots[1].set_ylabel("Pressure (Pa)")
    PT_plots[0].set_title("Temperature "+gas_label)
    PT_plots[1].set_title("Pressure "+gas_label)
    PT_plots[0].legend()
    PT_plots[1].legend()
    return

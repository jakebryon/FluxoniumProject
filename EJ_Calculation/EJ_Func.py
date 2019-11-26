### Function for calculating EJ and EL given a set of input voltages and 
### currents. A plot is generated

import numpy as np
import matplotlib.pyplot as plt


### Values for FDH03 C3
Current = np.array([0.03, 0.06, 0.1])

## EJ values
V200 = np.array([0.44, 0.68, 1.0])
V300 = np.array([0.32, 0.48, 0.67])
V400 = np.array([0.29, 0.40, 0.58])
V500 = np.array([0.26, 0.37, 0.51])
V600 = np.array([0.21, 0.30, 0.41])

V200_2 = np.array([0.55, 0.91, 1.47])
V300_2 = np.array([0.34, 0.55, 0.82])
V400_2 = np.array([0.28, 0.40, 0.60])
V500_2 = np.array([0.25, 0.35, 0.51])
V600_2 = np.array([0.21, 0.30, 0.43])


Volt_EJ = np.array([V200, V300, V400, V500, V600])
Volt_EJ_2 = np.array([V200_2, V300_2, V400_2, V500_2, V600_2])

widths = np.array([200.0,300.0,400.0,500.0,600.0])


## EL values
V40 = np.array([2.0, 3.75, 6.2])
V50 = np.array([2.56, 4.7, 7.8])
V60 = np.array([3.0, 5.6, 9.3])
V70 = np.array([3.4, 6.4, 10.5])

Volt_EL = np.array([V40, V50, V60, V70])
NumJunc = np.array([40, 50, 60, 70])



def EJ_Calc(Current, Voltages, widths, Voltages_2 = 0, EL = False):
    # Current is sourced from kethley and voltage is measured
    # this measurement will give us the resitances
    # Widths are widths of each junction
    # if doing EL calculation, working with number of undercuts, not widths

    # create empty list to fill in with resistances
    NumVolts = len(Voltages)
    Resistance_List = np.zeros(NumVolts)

    # Voltages_2 are the volttages from a second set of witnesses use for stats

    # find energies for different widths
    for i in range(NumVolts):
        EJ_Fit = np.polyfit(Voltages[i], Current, 1)
        EJ_Fit_Slope = EJ_Fit[0]
        Resistance_List[i] = 1/EJ_Fit_Slope

    Energies_List = 132.6/Resistance_List

    # if there is a second set of juncion measurements, find average
    if isinstance(Voltages_2, (list, tuple, np.ndarray)): 
        Resistance_List_2 = np.zeros(NumVolts)

        for i in range(NumVolts):
            EJ_Fit = np.polyfit(Voltages_2[i], Current, 1)
            EJ_Fit_Slope = EJ_Fit[0]
            Resistance_List_2[i] = 1/EJ_Fit_Slope

        Energies_List_2 = 132.6/Resistance_List_2
        Energies_List = (Energies_List + Energies_List_2)/2.0


    Energies_Fit_params = np.polyfit(widths, Energies_List, 1)
    Energies_Fit_widths = np.linspace(
                                    widths[0] - 50, 
                                    widths[-1] +50,
                                    100
                                    )

    Energies_Fit_Energies = (Energies_Fit_params[0]*Energies_Fit_widths
                            + Energies_Fit_params[1])

    plt.figure(1)
    plt.plot(widths, Energies_List, 'o')
    plt.plot(Energies_Fit_widths, Energies_Fit_Energies)
    plt.xlabel('Junction widths')
    plt.ylabel('Energy (GHz)')
    plt.show()

EJ_Calc(Current, Volt_EJ, widths, Voltages_2 = Volt_EJ_2)



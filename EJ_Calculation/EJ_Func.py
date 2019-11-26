### Function for calculating EJ and EL given a set of input voltages and 
### currents. A plot is generated

import numpy as np
import matplotlib.pyplot as plt


def EJ_Calc(
        Current, 
        Voltages, # array containing lists of measured voltages
        widths,  # width of junction or number of undercuts
        Voltages_2 = 0, # option for averaging with extra data
        FigNum = 1,  # number of figure created
        EL = False # tells if doing calcuation for EL not EJ
        ):
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

    plt.figure(FigNum)
    plt.plot(widths, Energies_List, 'o')
    plt.plot(Energies_Fit_widths, Energies_Fit_Energies)
    plt.xlabel('Junction widths')
    if EL:
        plt.xlabel('Number of undercuts')
    plt.ylabel('Energy (GHz)')
    plt.show()


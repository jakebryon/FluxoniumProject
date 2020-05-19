### Function for calculating EJ and EL given a set of input voltages and 
### currents. A plot is generated

import numpy as np
import matplotlib.pyplot as plt


# function for just calculating EJ
def EJ_Calc_basic(
    Current, # array of currents used for finding Voltages
    Voltages, # measured volttages for finding resistances
    ):
    # create empty list to fill in with resistances

    Res = (Voltages[1] - Voltages[0])/(Current[1] - Current[0])

    print('basic calculation of resistance in kOhms')
    print(Res)
    Energy = 132.6/Res
    return Energy


# function for calculating resitances
def Res_Calc(
    Current, # array of currents used for finding Voltages
    Voltages, # measured volttages for finding resistances
    ):
    # create empty list to fill in with resistances
    NumVolts = len(Voltages)
    Resistance_List = np.zeros(NumVolts)

    # find energies for different widths
    for i in range(NumVolts):
        EJ_Fit = np.polyfit(Voltages[i], Current, 1)
        EJ_Fit_Slope = EJ_Fit[0]
        Resistance_List[i] = 1/EJ_Fit_Slope

    return Resistance_List


# function for just calculating EJ
def EJ_Calc(
    Current, # array of currents used for finding Voltages
    Voltages, # measured volttages for finding resistances
    ):
    Resistance_List = Res_Calc(Current, Voltages) 
    
    Energies_List = 132.6/Resistance_List
    return Energies_List


# function for calculating and plotting Data
def EJ_Calc_Plot(
        Current, 
        Voltages, # array containing lists of measured voltages
        widths,  # width of junction or number of undercuts
        Voltages_2 = 0, # option for averaging with extra data
        FigNum = 1,  # number of figure created
        FigName = '', # name of plotted figure
        EL = False, # tells if doing calcuation for EL not EJ
        Desired_Energy = 0.0 # desired EJ/EL value, if this number is 
                        # given then numbers are given to achieve this
        ):

    NumVolts = len(Voltages)
    Energies_List = EJ_Calc(Current, Voltages)

    # setting extra range for EJ
    ext_range = 50

    # if there is a second set of juncion measurements, find average
    if isinstance(Voltages_2, np.ndarray): 
        print('test')
        
        Energies_List_2 = EJ_Calc(Current, Voltages_2)
        ### average together both lists of energies
        Energies_List = (Energies_List + Energies_List_2)/2.0
        ### setting fitting extra range for EL
        ext_range = 10


    # preform fit of energies
    Energies_Fit_params = np.polyfit(widths, Energies_List, 1)
    Energies_Fit_widths = np.linspace(
                                    widths[0] - ext_range, 
                                    widths[-1] + ext_range,
                                    100
                                    )

    Energies_Fit_Energies = (Energies_Fit_params[0]*Energies_Fit_widths
                            + Energies_Fit_params[1])

    # if given a desired Energy, find width or UnderCut
    if Desired_Energy != 0.0:

        Desired_Fit = np.array([
            Energies_Fit_params[0], 
            Energies_Fit_params[1] - Desired_Energy])
        
        Desired_Num = np.roots(Desired_Fit)
        Desired_var = 'width'
        if EL:
            Desired_var = 'Number of UnderCuts'
    
        print('Desired '+ Desired_var + ': ' + str(Desired_Num))


    # ploting data
    plt.figure(FigNum)
    plt.plot(widths, Energies_List, 'o')
    if not EL:
        plt.plot(Energies_Fit_widths, Energies_Fit_Energies)
    plt.xlabel('Junction widths')
    if EL:
        plt.xlabel('Number of undercuts')
    plt.ylabel('Energy (GHz)')
    plt.title(FigName)

def Energy_Stats(Energy_List):
    # gives array of avg, std, and string describing
    Energy_std = np.std(Energy_List)
    Energy_avg = np.mean(Energy_List)
    
    Energy_StatsStr = ('EJ = ' + 
        '{0:.2f}'.format(Energy_avg) + 
        u" \u00B1" + ' ' 
        '{0:.2f}'.format(Energy_std) 
        )
    return Energy_avg, Energy_std, Energy_StatsStr

def Energy_Stats_Plot(Energy_List, width):
    # plots error bar on an exisiting plot
    Energy_avg, Energy_std, Energy_StatsStr = Energy_Stats(Energy_List)

    plt.errorbar(width, Energy_avg, 
        yerr = Energy_std, 
        fmt = 'go',
        label = Energy_StatsStr
        )
    plt.legend()



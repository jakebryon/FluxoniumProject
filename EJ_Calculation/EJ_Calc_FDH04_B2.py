### calculating EJ and EL using functions in EJ_Func.py
### FDH04_B2
### This device was a dose test device

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from EJ_Func import *

### Values for FDH03 B2
Current = np.array([0.10, 0.30, 0.50])

### finding Energy of EJ junction, width = 225.0
CurrentEJ = np.array([0.10, 0.50])
EJ_Volt = np.array([2.17, 10.3])
EJ_Volt_2 = EJ_Volt
EJs_Volt = np.array([EJ_Volt, EJ_Volt_2])
EJs_Energy_List = EJ_Calc(CurrentEJ, EJs_Volt)
EJs_width = 225.0

print('EJ Res')
EJ_Calc_basic(Current, EJ_Volt)

EJs_Energy_avg, EJs_Energy_std, EJs_StatsStr = Energy_Stats(EJs_Energy_List)


# finding Energy of EL junction, UnderCuts = 100
CurrentEL = np.array([0.1, 0.5, 1.0])
EL_1_Volt = np.array([26.0, 126.0, 251.0])
EL_2_Volt = np.array([23.2, 114.0, 228.0])
ELs_Volt = np.array([EL_1_Volt, EL_2_Volt])
ELs_Res_List = Res_Calc(Current, ELs_Volt)
ELs_Energy_List = EJ_Calc(Current, ELs_Volt)
ELs_UnderCut = 100.0

print('EL res')
print(ELs_Res_List)


ELs_Energy_avg, ELs_Energy_std, ELs_StatsStr = Energy_Stats(ELs_Energy_List)


### Fluxonium test junction resistance
FLX_Current = np.array([0.1, 0.5, 1.0])
FLX_Volt = np.array([2.0, 9.84, 19.6])
FLX_Volt_2 = np.array([5.3, 15.8, 26.2])
EJ_Calc_basic(FLX_Current, FLX_Volt)
EJ_Calc_basic(FLX_Current, FLX_Volt_2)



########################## Wittness EJ values
V100 = np.array([5.3, 15.5, 25.7])
V200 = np.array([1.60, 4.60, 7.60])
V300 = np.array([1.23, 3.52, 5.80])
V400 = np.array([0.93, 2.60, 4.30])
V500 = np.array([0.74, 2.05, 3.35])

V100_2 = np.array([2.56, 7.53, 12.4])
V200_2 = np.array([1.40, 4.05, 6.70])
V300_2 = np.array([1.12, 3.23, 5.32])
V400_2 = np.array([0.86, 2.44, 4.00])
V500_2 = np.array([0.695, 1.95, 3.20])

Volt_EJ = np.array([V100, V200, V300, V400, V500])
Volt_EJ_2 = np.array([V100_2, V200_2, V300_2, V400_2, V500_2])
widths = np.array([100.0, 200.0, 300.0, 400.0, 500.0])

### removing the 100 length
Volt_EJ = np.array([V200, V300, V400, V500])
Volt_EJ_2 = np.array([V200_2, V300_2, V400_2, V500_2])
widths = np.array([200.0, 300.0, 400.0, 500.0])



EJ_Res = Res_Calc(Current, Volt_EJ)
EJ_Res_2 = Res_Calc(Current, Volt_EJ_2)
EJ_Eng = EJ_Calc(Current, Volt_EJ)
EJ_Eng_2 = EJ_Calc(Current, Volt_EJ_2)

EJ_avg1 = 132.6/((EJ_Res + EJ_Res_2)/2)
EJ_avg2 = (EJ_Eng + EJ_Eng_2)/2
print('EJ energy list')
print(EJ_avg1)
print(EJ_avg2)

### EJ_Calc(Current, Volt_EJ)

## Wittness EL values
V10 = np.array([2.15, 6.27, 10.3])
V20 = np.array([4.30, 12.8, 21.1])
V30 = np.array([6.50, 19.2, 31.8])
V40 = np.array([8.80, 26.2, 43.3])
V50 = np.array([10.7, 31.8, 52.6])

V10_2 = np.array([2.03, 5.96, 9.86])
V20_2 = np.array([4.00, 11.8, 19.6])
V30_2 = np.array([6.09, 18.1, 30.1])
V40_2 = np.array([8.10, 24.2, 40.2])
V50_2 = np.array([10.7, 31.8, 52.6])

# adding in points for EL
Volt_EL = np.array([V10, V20, V30, V40, V50])
Volt_EL_2 = np.array([V10_2, V20_2, V30_2, V40_2, V50_2])
NumJunc = np.array([10, 20, 30, 40, 50])

EJ_Calc(Current, Volt_EL)

EL_Eng = EJ_Calc(Current, Volt_EL)
EL_Eng_2 = EJ_Calc(Current, Volt_EL_2)
EL_Energy = (EL_Eng + EL_Eng_2)/2



############################################## run calculation
# printing basic Stats
print(EJs_StatsStr)
print(ELs_StatsStr)
####################################### plotting EJ for FDH03
EJ_Calc_Plot(
    Current, Volt_EJ, widths, 
    Voltages_2 = Volt_EJ_2,
    FigNum = 1,
    FigName = 'FDH04_B2 Wittness Junctions: EJ',
    Desired_Energy = 8.35
    )
# EJ_Calc_Plot(
#     Current, Volt_EJ_2, widths, 
# #    Voltages_2 = Volt_EJ_2,
#     FigNum = 1,
#     FigName = 'FDH04_B2 Wittness Junctions: EJ, part 2',
#     Desired_Energy = 8.35
#     )
Energy_Stats_Plot(EJs_Energy_List, EJs_width)

###################################### plotting EL for FDH03
EJ_Calc_Plot(
    Current, Volt_EL, NumJunc,
    Voltages_2 = Volt_EL_2, 
    FigNum = 2, 
    FigName = 'FDH04_B2 Wittness Junctions EL',
    EL = True,
    Desired_Energy = 0.45
    )
Energy_Stats_Plot(ELs_Energy_List, ELs_UnderCut)

plt.show()



######################### testing 1/x fit for EL

x_data = np.array([10, 20, 30, 40, 50, 100])
y_data = np.array([
    6.64081207, 3.27849686, 2.1527556,  1.59484183, 1.26584949, 0.25788128
    ])

# x_data = np.array([10, 20, 30, 40, 50, 100])
# y_data = np.append(EL_Energy, ELs_Energy_List[0])

print(y_data)

# print(ELs_Energy_List)

def func(x, a, c):
    return (a/(x)) + c

popt, pcov = curve_fit(func, x_data, y_data)

x_data_lin = np.linspace(8, 105, 100)

plt.figure(10)
plt.plot(x_data_lin, func(x_data_lin, *popt), 'r-')
plt.plot(x_data, y_data, 'o')

plt.title('FDH04 EL fit')
plt.xlabel('number of undercuts')
plt.ylabel('Energy (GHz)')

plt.show()

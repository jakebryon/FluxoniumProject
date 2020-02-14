### calculating EJ and EL using functions in EJ_Func.py
### FDH04_C2
### This device was a dose test device

import numpy as np
import matplotlib.pyplot as plt
from EJ_Func import *

### Values for FDH03 C2
Current = np.array([0.05, 0.10, 0.15])

# finding Energy of EJ junction, width = 225.0
CurrentEJ = np.array([0.05, 0.10, 0.15])
EJ_Volt = np.array([0.54, 0.95, 1.36])
EJ_Volt_2 = np.array([0.54, 0.95, 1.36])
EJs_Volt = np.array([EJ_Volt, EJ_Volt_2])
EJs_Energy_List = EJ_Calc(CurrentEJ, EJs_Volt)
EJs_width = 225.0

### EJ_Calc_basic(Current, EJ_Volt)

EJs_Energy_avg, EJs_Energy_std, EJs_StatsStr = Energy_Stats(EJs_Energy_List)


# finding Energy of EL junction, UnderCuts = 100
EL_1_Volt = np.array([9.80,18.6,27.4])
EL_2_Volt = np.array([10.5,19.9,23.3])
ELs_Volt = np.array([EL_1_Volt, EL_2_Volt])
ELs_Energy_List = EJ_Calc(Current, ELs_Volt)
ELs_UnderCut = 100.0

ELs_Energy_avg, ELs_Energy_std, ELs_StatsStr = Energy_Stats(ELs_Energy_List)


########################## Wittness EJ values
V200 = np.array([0.62, 1.14, 1.65])
V225 = np.array([0.54, 0.95, 1.36])
V300 = np.array([0.46, 0.82, 1.19])
V500 = np.array([0.35, 0.63, 0.90])

Volt_EJ = np.array([V200, V225, V300, V500])
widths = np.array([200.0,225.0,300.0,500.0])

### EJ_Calc(Current, Volt_EJ)

## Wittness EL values
V10 = np.array([0.95, 1.73, 2.54])
V20 = np.array([1.93, 3.59, 5.28])
V30 = np.array([2.85, 5.42, 8.03])
V40 = np.array([3.94, 7.51, 11.1])
V50 = np.array([5.20, 9.55, 14.0])

V10_2 = np.array([1.00, 1.85, 2.70])
V20_2 = np.array([2.13, 3.89, 5.69])
V30_2 = np.array([3.10, 5.91, 8.74])
V40_2 = np.array([4.30, 8.08, 11.9])
V50_2 = np.array([5.00, 9.58, 14.2])

# adding in points for EL
Volt_EL = np.array([V10, V20, V30, V40, V50, EL_1_Volt])
Volt_EL_2 = np.array([V10_2, V20_2, V30_2, V40_2, V50_2, EL_2_Volt])
NumJunc = np.array([10, 20, 30, 40, 50, ELs_UnderCut])

EJ_Calc(Current, Volt_EL)


############################################## run calculation
# printing basic Stats
print(EJs_StatsStr)
print(ELs_StatsStr)
####################################### plotting EJ for FDH03
EJ_Calc_Plot(
    Current, Volt_EJ, widths, 
    FigNum = 1,
    FigName = 'FDH04_C2 Wittness Junctions: EJ',
    Desired_Energy = 8.35
    )
Energy_Stats_Plot(EJs_Energy_List, EJs_width)

###################################### plotting EL for FDH03
EJ_Calc_Plot(
    Current, Volt_EL, NumJunc,
    Voltages_2 = Volt_EL_2, 
    FigNum = 2, 
    FigName = 'FDH04_C2 Wittness Junctions EL',
    EL = True,
    Desired_Energy = 0.45
    )
Energy_Stats_Plot(ELs_Energy_List, ELs_UnderCut)

plt.show()

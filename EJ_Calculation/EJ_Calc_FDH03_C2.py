### calculating EJ and EL using functions in EJ_Func.py
### FDH03_C2
### When making the test junctions for this device, I forgot to change
### the EJ and EL wittness to match the actualy ones
## EL should be 135
## EJ should be 120

import numpy as np
import matplotlib.pyplot as plt
from EJ_Func import *

### Values for FDH03 C2
Current = np.array([0.05, 0.10, 0.15])

# finding Energy of EJ junction, width = 350.0
CurrentEJ = np.array([0.05, 0.10, 0.15])
EJ_1_Volt = np.array([0.54, 0.95, 1.36])
EJ_2_Volt = np.array([0.60, 1.12, 1.63])
EJ_3_Volt = np.array([0.50, 0.91, 1.30])
EJ_4_Volt = np.array([0.52, 0.97, 1.40])
EJs_Volt = np.array([EJ_1_Volt, EJ_2_Volt, EJ_3_Volt, EJ_4_Volt])
EJs_Energy_List = EJ_Calc(CurrentEJ, EJs_Volt)
EJs_width = 350.0

EJs_Energy_avg, EJs_Energy_std, EJs_StatsStr = Energy_Stats(EJs_Energy_List)


# finding Energy of EL junction, UnderCuts = 120
EL_1_Volt = np.array([12.0,22.7, 33.5])
EL_2_Volt = np.array([11.4,22.0, 32.6])
EL_3_Volt = np.array([12.0,23.2, 34.4])
EL_4_Volt = np.array([11.5,22.4, 33.4])
ELs_Volt = np.array([EL_1_Volt, EL_2_Volt, EL_3_Volt, EL_4_Volt])
ELs_Energy_List = EJ_Calc(Current, ELs_Volt)
ELs_UnderCut = 120.0

ELs_Energy_avg, ELs_Energy_std, ELs_StatsStr = Energy_Stats(ELs_Energy_List)


########################## Wittness EJ values
V200 = np.array([1.0, 1.88, 2.74])
V300 = np.array([0.61, 1.12, 1.61])
V400 = np.array([0.50, 0.88, 1.27])
V500 = np.array([0.33, 0.57, 0.79])
V600 = np.array([0.26, 0.45, 0.62])

V200_2 = np.array([1.37, 2.62, 3.85])
V300_2 = np.array([0.80, 1.50, 2.19])
V400_2 = np.array([0.52, 0.95, 1.38])
V500_2 = np.array([0.35, 0.61, 0.90])
V600_2 = np.array([0.30, 0.46, 0.64])

Volt_EJ = np.array([V200, V300, V400, V500, V600])
Volt_EJ_2 = np.array([V200_2, V300_2, V400_2, V500_2, V600_2])
widths = np.array([200.0,300.0,400.0,500.0,600.0])

## Wittness EL values
V40 = np.array([3.5, 6.63, 9.8])
V50 = np.array([4.3, 8.1, 12.1])
V60 = np.array([5.2, 10.0, 14.7])
V70 = np.array([6.0, 11.5, 17.1])

V40_2 = np.array([3.5, 6.7, 9.9])
V50_2 = np.array([4.1, 7.9, 11.8])
V60_2 = np.array([5.27, 10.2, 15.2])
V70_2 = np.array([6.0, 11.7, 17.5])

# adding in points for EL
Volt_EL = np.array([V40, V50, V60, V70, EL_1_Volt])
Volt_EL_2 = np.array([V40_2, V50_2, V60_2, V70_2, EL_2_Volt])
NumJunc = np.array([40, 50, 60, 70, ELs_UnderCut])


############################################## run calculation
# printing basic Stats
print(EJs_StatsStr)
print(ELs_StatsStr)
####################################### plotting EJ for FDH03
EJ_Calc_Plot(
    Current, Volt_EJ, widths, 
    Voltages_2 = Volt_EJ_2, 
    FigNum = 1,
    FigName = 'FDH03_C3 Wittness Junctions: EJ',
    Desired_Energy = 8.35
    )
Energy_Stats_Plot(EJs_Energy_List, EJs_width)

###################################### plotting EL for FDH03
EJ_Calc_Plot(
    Current, Volt_EL, NumJunc,
    Voltages_2 = Volt_EL_2, 
    FigNum = 2, 
    FigName = 'FDH03_C3 Wittness Junctions EL',
    EL = True,
    Desired_Energy = 0.45
    )
Energy_Stats_Plot(ELs_Energy_List, ELs_UnderCut)

plt.show()

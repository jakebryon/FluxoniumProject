### calculating EJ and EL using functions in EJ_Func.py
### FDH03_C3

import numpy as np
import matplotlib.pyplot as plt
from EJ_Func import *

### Values for FDH03 C3
Current = np.array([0.03, 0.06, 0.1])

## Wittness EJ values
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

# finding Energy of EJ junction, width = 350
EJ_1_Volt = np.array([0.32, 0.46, 0.64])
EJ_2_Volt = np.array([0.30, 0.44, 0.64])
EJ_3_Volt = np.array([0.34, 0.53, 0.79])
EJ_4_Volt = np.array([0.30, 0.45, 0.67])
EJs_Volt = np.array([EJ_1_Volt, EJ_2_Volt, EJ_3_Volt, EJ_4_Volt])
EJs_Energy_List = EJ_Calc(Current, EJs_Volt)
EJs_width = 350.0

EJs_Energy_avg, EJs_Energy_std, EJs_StatsStr = Energy_Stats(EJs_Energy_List)

## Wittness EL values
#V40 = np.array([2.0, 3.75, 6.2])
V50 = np.array([2.56, 4.7, 7.8])
V60 = np.array([3.0, 5.6, 9.3])
V70 = np.array([3.4, 6.4, 10.5])

#V40_2 = np.array([1.9, 3.6, 6.0])
V50_2 = np.array([2.5, 4.7, 7.8])
V60_2 = np.array([2.9, 5.6, 9.4])
V70_2 = np.array([3.4, 6.5, 11.0])

# adding in points for EL
EL_1_Volt = np.array([7.0, 12.5, 20.3])
EL_2_Volt = np.array([7.1, 13.4, 22.4])

#Volt_EL = np.array([V40, V50, V60, V70, EL_1_Volt])
#Volt_EL_2 = np.array([V40_2, V50_2, V60_2, V70_2, EL_2_Volt])
#NumJunc = np.array([40, 50, 60, 70, 122])

Volt_EL = np.array([V50, V60, V70, EL_1_Volt])
Volt_EL_2 = np.array([V50_2, V60_2, V70_2, EL_2_Volt])
NumJunc = np.array([50, 60, 70, 122])



# finding Energy of EL junction, UnderCuts = 120
EL_3_Volt = np.array([6.0, 11.4, 18.8])
EL_4_Volt = np.array([6.5, 12.5, 21.1])
ELs_Volt = np.array([EL_1_Volt, EL_2_Volt, EL_3_Volt, EL_4_Volt])
ELs_Energy_List = EJ_Calc(Current, ELs_Volt)
ELs_UnderCut = 120.0

ELs_Energy_avg, ELs_Energy_std, ELs_StatsStr = Energy_Stats(ELs_Energy_List)

#### run calculation
# printing basic Stats
print(EJs_StatsStr)
print(ELs_StatsStr)
# plotting EJ for FDH03
EJ_Calc_Plot(
    Current, Volt_EJ, widths, 
    Voltages_2 = Volt_EJ_2, 
    FigNum = 1,
    FigName = 'FDH03_C3 Wittness Junctions: EJ',
    Desired_Energy = 8.35
    )
Energy_Stats_Plot(EJs_Energy_List, EJs_width)

# plotting EL for FDH03
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

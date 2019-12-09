### calculating EJ and EL using functions in EJ_Func.py
### FDH02
import numpy as np
import matplotlib.pyplot as plt
from EJ_Func import *


### Values for FDH02 B2
Current = np.array([0.03, 0.06, 0.1])
### Neglecting values for 200 and 250 because they seem off
### by off I mean non-linear with the others
#V200 = np.array([3.8, 7.4, 12.4])
V250 = np.array([3.8, 7.6, 12.6])
V300 = np.array([0.85, 1.6, 2.6])
V350 = np.array([0.55, 1.0, 1.6])
V400 = np.array([0.40, 0.75, 1.2])
V450 = np.array([0.35, 0.60, 0.95])
V500 = np.array([0.31, 0.55, 0.87])
V550 = np.array([0.29, 0.5, 0.8])
V600 = np.array([0.27, 0.47, 0.74])

Volt_EJ = np.array([V250, V300, V350, V400, V450, V500, V550, V600])

# widths of junctions
widths = np.array([250,300,350,400,450,500,550,600])

EJs_Energy_List = EJ_Calc(Current, Volt_EJ)

EJs_Energy_avg, EJs_Energy_std, EJs_StatsStr = Energy_Stats(EJs_Energy_List)

############################################## run calculation
# printing basic Stats
print(EJs_StatsStr)
####################################### plotting EJ for FDH02_B2
EJ_Calc_Plot(
    Current, Volt_EJ, widths, 
    FigNum = 1,
    FigName = 'FDH02_B2 Wittness Junctions: EJ',
    )

plt.show()

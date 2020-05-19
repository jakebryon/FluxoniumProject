### Code for calculating EJ given junction resistances

import numpy as np
import matplotlib.pyplot as plt


# ### Values for FDH02 B2
# Current = np.array([0.03, 0.06, 0.1])
# ### Neglecting values for 200 and 250 because they seem off
# ### by off I mean non-linear with the others
# #V200 = np.array([3.8, 7.4, 12.4])
# V250 = np.array([3.8, 7.6, 12.6])
# V300 = np.array([0.85, 1.6, 2.6])
# V350 = np.array([0.55, 1.0, 1.6])
# V400 = np.array([0.40, 0.75, 1.2])
# V450 = np.array([0.35, 0.60, 0.95])
# V500 = np.array([0.31, 0.55, 0.87])
# V550 = np.array([0.29, 0.5, 0.8])
# V600 = np.array([0.27, 0.47, 0.74])
# 
# Volt = np.array([V250, V300, V350, V400, V450, V500, V550, V600])
# 
# # widths are given by x*2000 nm, where x is width of junction and is
# # the value in the following array, the array below represents 
# # widths of junctions
# widths = np.array([250,300,350,400,450,500,550,600])
# Resistances = np.zeros(len(Volt))
# Intercepts = np.zeros(len(Volt))
#### ######################################################################
#### 
#### 
#### for i in range(0, len(Resistances)):
####     slope = np.polyfit(Volt[i], Current, 1)[0]
####     Resistances[i] = 1/slope
#### 
#### #Energies = 132.6*Resistances
#### Energies = 132.6/Resistances
#### 
#### # using Energy and area, find linear approx;
#### Energy_params = np.polyfit(widths, Energies, 1)
#### widths_fit = np.linspace(150, 660, 300)
#### Energies_fit = widths_fit*Energy_params[0] + Energy_params[1]
#### 
#### # finding Energy of FDH02 EJ
#### FDH02_A = 190.0
#### FDH02_E = FDH02_A*Energy_params[0] + Energy_params[1]
#### print('FDH02 width: ', FDH02_A, 'FDH02 EJ: ', FDH02_E, 'GHz')
#### 
#### # Note: On Witness Junctions, There was no copied over EL, no EL value!
#### # however, there was an old EL chain from FDH01
#### FDH01_V = np.array([3.5, 7.0, 11.5])
#### slope = np.polyfit(FDH01_V, Current, 1)[0]
#### FDH01_R = 1/slope
#### FDH01_EL = 132.6/FDH01_R
#### print('FDH01 EL ', FDH01_EL)
#### 
#### 
#### 
#### plt.figure(1)
#### plt.plot(widths, Energies, '-o')
#### plt.plot(widths_fit, Energies_fit, '-')
#### plt.plot(FDH02_A, FDH02_E, '*')
#### plt.xlabel('Junction width (nm)')
#### plt.ylabel('Energy (GHz)')
#### plt.title('FDH02 Witness Junction Energies')
#### plt.legend(['Witness Junctions','Fit','FDH02 EJ from fit'])
#### 
#### plt.show()
#### 
#### 
##########################################################################

### Values for FDH03 C3
Current = np.array([0.03, 0.06, 0.1])

V200 = np.array([0.44, 0.68, 1.0])
V300 = np.array([0.32, 0.48, 0.67])
V400 = np.array([0.29, 0.40, 0.58])
V500 = np.array([0.26, 0.37, 0.51])
V600 = np.array([0.21, 0.30, 0.41])

Volt_EJ = np.array([V200, V300, V400, V500, V600])
widths = np.array([200,300,400,500,600])


V40 = np.array([2.0, 3.75, 6.2])
V50 = np.array([2.56, 4.7, 7.8])
V60 = np.array([3.0, 5.6, 9.3])
V70 = np.array([3.4, 6.4, 10.5])

Volt_EL = np.array([V40, V50, V60, V70])
NumJunc = np.array([40, 50, 60, 70])


EJ_1 = np.array([0.32, 0.46, 0.64])
EJ_2 = np.array([0.30, 0.44, 0.64])
EJs = np.array([EJ_1, EJ_2])

EL_1 = np.array([7.0, 12.5, 20.3])
EL_2 = np.array([7.1, 13.4, 22.4])
ELs = np.array([EL_1, EL_2])


# widths are given by x*2000 nm, where x is width of junction and is
# the value in the following array, the array below represents 
# widths of junctions

Resistances_EJ = np.zeros(len(Volt_EJ))
Resistances_EL = np.zeros(len(Volt_EL))
Resistances_EJs = np.zeros(len(EJs))
Resistances_ELs = np.zeros(len(ELs))

######################################################################


for i in range(0, len(Resistances_EJ)):
    slope_EJ = np.polyfit(Volt_EJ[i], Current, 1)[0]
    Resistances_EJ[i] = 1/slope_EJ
print(Resistances_EJ)

for i in range(0, len(Resistances_EL)):
    slope_EL = np.polyfit(Volt_EL[i], Current, 1)[0]
    Resistances_EL[i] = 1/slope_EL

for i in range(0, len(Resistances_EJs)):
    slope_EJs = np.polyfit(EJs[i], Current, 1)[0]
    Resistances_EJs[i] = 1/slope_EJs

for i in range(0, len(Resistances_ELs)):
    slope_ELs = np.polyfit(ELs[i], Current, 1)[0]
    Resistances_ELs[i] = 1/slope_ELs



Energies_EJ = 132.6/Resistances_EJ
Energies_EL = 132.6/Resistances_EL
Energy_EJ = 132.6/Resistances_EJs
Energy_EL = 132.6/Resistances_ELs

EJ_Energy = np.average(Energy_EJ)
EL_Energy = np.average(Energy_EL)


# using Energy and area, find linear approx;
Energy_EJ_params = np.polyfit(widths, Energies_EJ, 1)
widths_fit = np.linspace(150, 650, 300)
Energies_EJ_fit = widths_fit*Energy_EJ_params[0] + Energy_EJ_params[1]

Energy_EL_params = np.polyfit(NumJunc, Energies_EL, 1)
NumJunc_fit = np.linspace(30, 80, 300)
Energies_EL_fit = NumJunc_fit*Energy_EL_params[0] + Energy_EL_params[1]

# printing out EJ and EC
print('EJ = ', EJ_Energy)
print('EL = ', EL_Energy)


plt.figure(1)
plt.suptitle('FDH03 Junction Energies')

ax1 = plt.subplot(1,2,1)
ax1.plot(widths, Energies_EJ, '-o')
ax1.plot(widths_fit, Energies_EJ_fit, '-')
ax1.set_xlabel('Junction width (nm)')
ax1.set_ylabel('Energy (GHz)')
ax1.set_title('Junction EJ')

ax2 = plt.subplot(1,2,2)
ax2.plot(NumJunc, Energies_EL, '-o')
ax2.plot(NumJunc_fit, Energies_EL_fit, '-')
ax2.set_xlabel('Number of Undercuts')
ax2.set_ylabel('Energy (GHz)')
ax2.set_title('Junction EL')

plt.show()


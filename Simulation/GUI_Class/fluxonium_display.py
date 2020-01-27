## 
# fluxonium_display contains the Display class which takes in a fluxonium
# object from scqubits and handles the plotting and display of 
# the GUI interface
# 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button
import os
from datetime import date
import csv
from analysis import *
import scqubits as qubit

class Display:
    ## initialize object by creating a default fluxonium object with
    # note: all energies are in GHz
    # input parameters
    # @param self the object pointer
    # @param EJ:float josphson junction energy
    # @param EC:float charging energy 
    # @param EL:float inductive energy
    # @param cav_freq: float cavity frequency
    # @param flux:float external magnetic flux
    # @param Num_levels:int number of levels to be in the fluxonium object
    # @parm fluxonium: qubit object of the simulated fluxonium object 

    def __init__(self, EJ, EC, EL, cav_freq, flux, Num_levels):
        ### EJ: float of JJ energy
        self.EJ = EJ
        ### EC: float of charging energy
        self.EC = EC
        ### EL: float of inductive energy
        self.EL = EL
        ### cav_freq: float of cavity frequency
        self.cav_freq = cav_freq
        ### flux: float of external flux
        self.flux = flux
        ### Num_levels: int of number of levels in model
        self.Num_levels = Num_levels

        self.fig, self._subplots_temp = ( 
            plt.subplots(2,2, figsize = [14.0, 5.0]) )
        #  turn format of subplots into array for simplicity 
        self.subplots = np.array([
           self._subplots_temp[0,0],
           self._subplots_temp[0,1],
           self._subplots_temp[1,0],
           self._subplots_temp[1,1]
           ])

        ## fluxonium: qubit object of the fluxonium
        self.fluxonium = qubit.Fluxonium(
            EJ = self.EJ,
            EC = self.EC,
            EL = self.EL,
            flux = self.flux,
            cutoff = 100
            )


    ## method to place down all the buttons/text entry boxes
    #
    def place_buttons(self):

        # adjust all subplots to make room for buttons
        self.fig.subplots_adjust(
            left = 0.18, 
            bottom = 0.1, 
            wspace = 0.2,
            hspace = 0.3,
            )
        self.fig.suptitle('Alira Laser Data Plotting')

        # placing EJ text entry
        axtext_EJ = plt.axes([0.05, 0.80, 0.06, 0.05])
        self.text_EJ = TextBox(axtext_EJ, 'EJ:', initial = str(self.EJ))
        self.text_EJ.on_submit(self.EJ_change)

        # placing EC text entry
        axtext_EC = plt.axes([0.05, 0.70, 0.06, 0.05])
        self.text_EC = TextBox(axtext_EC, 'EC:', initial = str(self.EC))
        self.text_EC.on_submit(self.EC_change)

        # placing EL text entry
        axtext_EL = plt.axes([0.05, 0.60, 0.06, 0.05])
        self.text_EL = TextBox(axtext_EL, 'EL:', initial = str(self.EL))
        self.text_EL.on_submit(self.EL_change)

        # placing cavity frequencu text entry
        axtext_cav = plt.axes([0.05, 0.50, 0.06, 0.05])
        self.text_cav = TextBox(axtext_cav, 'cavity:', 
            initial = str(self.cav_freq))
        self.text_cav.on_submit(self.cavity_change)

        # placing phi text entry
        axtext_phi = plt.axes([0.05, 0.40, 0.06, 0.05])
        self.text_phi = TextBox(axtext_phi, 'phi_ext', initial = str(self.flux))
        self.text_phi.on_submit(self.phi_change)


    ## method to plot the energy levels vs external flux
    #
    def flux_sweep_plot(self, plot_num = 0):
        # if method is called set to true
        self.flux_sweep = True
        # keep track of plot number
        self.flux_sweep_num = plot_num

        self.subplots[plot_num].cla()
        phi_ext = np.linspace(0.0, 1.0, 100)
        
        energy_spec = self.fluxonium.get_spectrum_vs_paramvals(
            'flux',
            phi_ext,
            evals_count = self.Num_levels
            )

        for i in range(0, self.Num_levels):
            self.subplots[plot_num].plot(
                energy_spec.param_vals,
                energy_spec.energy_table[:,i] - energy_spec.energy_table[:,0],
                '-'
                )

        self.subplots[plot_num].set_xlabel(r'$\varphi_{ext}/2\pi$')
        self.subplots[plot_num].set_ylabel('E - E(0) (GHz)')
#        self.subplots[plot_num].legend(legend)
        self.subplots[plot_num].set_title(
            'Energy Levels vs. ' + r'$\varphi_{ext}/2\pi$')

    ## method to change the value for EJ, used in text box input
    #
    def EJ_change(self, text):
       self.EJ = float(text)
       self.fluxonium.EJ = self.EJ
       self.update()


    ## method to change the value for EC
    #
    def EC_change(self, text):
       self.EC = float(text)
       self.fluxonium.EC = self.EC
       self.update()


    ## method to change the value for EL
    #
    def EL_change(self, text):
       self.EL = float(text)
       self.fluxonium.EC = self.EL
       self.update()


    ## method to change the value for EL
    #
    def cavity_change(self, text):
       self.cav_freq = float(text)
       self.update()


    ## method to change the value for EL
    #
    def phi_change(self, text):
       self.flux = float(text)
       self.fluxonium.flux = self.flux
       self.update()


    ## method to update all running plots
    #
    def update(self):
        if self.flux_sweep == True:
            self.flux_sweep_plot(plot_num = self.flux_sweep_num)


    ## method to show plots
    #
    def show(self):
        plt.show()


#########################################################
### testing with a quick driver

display = Display(
    EJ = 8.34,
    EC = 1.15,
    EL = 0.45,
    cav_freq = 6.0,
    flux = 0.75, 
    Num_levels = 5
    )

display.place_buttons()
display.flux_sweep_plot()
display.show()


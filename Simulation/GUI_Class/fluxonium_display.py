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
    # @pram Num_sum: int number of levels used in sum calculations
    # @parm fluxonium: qubit object of the simulated fluxonium object 

    def __init__(
            self, 
            EJ, 
            EC, 
            EL, 
            cav_freq, 
            flux, 
            Num_levels, 
            Num_sum = 10,
            Num_plots = 4,
            flux_start = 0.0,
            flux_stop = 1.0
            ):

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
        ### Num_sum: int of the number of levels to be used in calculations
        self.Num_sum = Num_sum
        ### Num_plots: int of the number plots on figure
        self.Num_plots = Num_plots
        ### flux_start is the starting flux for flux sweep
        self.flux_start = flux_start
        ### flux_stop is the stopping flux for flux sweep
        self.flux_stop = flux_stop


        # creating subplots given number requested
        # for 1-3 plots, just lay out in a row
        if self.Num_plots == 1:
            self.fig, self.subplots = ( 
                plt.subplots(1,1, figsize = [14.0, 5.0]) )
            ## subplots: list containing all subplots of the figure
            self.subplots = [self.subplots]

        elif self.Num_plots == 2:
            self.fig, self.subplots = ( 
                plt.subplots(1,2, figsize = [14.0, 5.0]) )

        elif self.Num_plots == 3:
            self.fig, self.subplots = ( 
                plt.subplots(1,3, figsize = [14.0, 5.0]) )

        elif self.Num_plots == 4: 
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

        ## set all booleans to false
        self.flux_sweep = False
        self.g_mat_bool = False
        self.wavefunciton_bool = False
        self.Raman_bool = False


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
        self.fig.suptitle('Fluxonium Interactive Plotting')

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


    ## method to plot the matrix elements vs external flux
    #
    def Flux_sweep_plot(self, plot_num = 0):
        # if method is called set to true
        self.flux_sweep = True
        # keep track of plot number
        self.flux_sweep_num = plot_num

        self.subplots[plot_num].cla()
        phi_ext = np.linspace(self.flux_start, self.flux_stop, 100)
        
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

        # plot cavity frequency

        self.subplots[plot_num].set_xlabel(r'$\varphi_{ext}/2\pi$')
        self.subplots[plot_num].set_ylabel('E - E(0) (GHz)')
#        self.subplots[plot_num].legend(legend)
        self.subplots[plot_num].set_title(
            'Energy Levels vs. ' + r'$\varphi_{ext}/2\pi$')


    ## method to plot the energy levels vs external flux
    # @param coupling_level level for which elements are calculted with
    # 
    def g_mat_plot(self, 
            plot_num = 0, 
            coupling_level = 0, 
            beta_phi = 0.2,
            fontsize = 10
            ):
        # if method is called set to true
        self.g_mat_bool = True
        # keep track of plot number
        self.g_mat_num = plot_num

        self.subplots[plot_num].cla()
        phi_ext = np.linspace(0.0, 1.0, 100)
        
        g_mat_array = self.fluxonium.get_matelements_vs_paramvals(
            operator = 'n_operator',
            param_name = 'flux',
            param_vals = phi_ext,
            evals_count = self.Num_levels
            )

        # g_const is the coupling constant factor
        g_const = self.cav_freq * beta_phi * (2.0/np.pi)
        g_const = 1 # just finding matrix elements

        g_val_array = np.zeros([len(phi_ext), self.Num_levels])
        for i in range(0, self.Num_levels):
            g_val_array[:, i] = np.abs(
                g_mat_array.matrixelem_table[:,coupling_level,i] )
            # zero out the self couplings
            if i == coupling_level: 
                g_val_array[:,i] = np.zeros(len(phi_ext))

        g_val_array = g_val_array*g_const

        for i in range(0, self.Num_levels):
            self.subplots[plot_num].plot(
                g_mat_array.param_vals,
                g_val_array[:,i], # *1e3 MHz units
                '-',
                label = 'level: '+str(i)
                )

        self.subplots[plot_num].set_xlabel(r'$\varphi_{ext}/2\pi$',
            fontsize=fontsize)
        ### defining y axis name
        yaxis = r'$| \langle$'+str(coupling_level)+r'$|\hat{n}| i \rangle |$'
        self.subplots[plot_num].set_ylabel(yaxis, 
            fontsize=fontsize)
        self.subplots[plot_num].set_title(
            'coupling level ' + str(coupling_level), 
            fontsize=fontsize)
        self.subplots[plot_num].legend()


    ## method to plot wave functions on the hamiltonian
    # @param plot_num int for which plot to be displayed on
    # 
    def wavefunction_plot(self, 
            plot_num = 0, 
            phi_limit = 3
            ):
        # if method is called set to true
        self.wavefunction_bool = True
        # keep track of plot number
        self.wavefunction_num = plot_num

        self.subplots[plot_num].cla()

        eigen_vals, eigen_vec = (
            self.fluxonium.eigensys(evals_count = self.Num_levels) )
        esys = (
            self.fluxonium.eigensys(evals_count = self.Num_levels) )

        wavefunc_array = np.array([0,0,0,0,0,0], dtype = object)
        wavefunc_array[0] = self.fluxonium.wavefunction(esys)

        for i in range(0, self.Num_levels):
            wavefunc_array[i] = self.fluxonium.wavefunction( 
                esys = esys,
                which = i,
                phi_range = (-phi_limit*np.pi, phi_limit*np.pi)
#                100
                )
        # plot wavefunctions
        for i in range(0, self.Num_levels):
            wavefunc = wavefunc_array[i]
            energy_scaled = eigen_vals[i] - eigen_vals[0]
            self.subplots[plot_num].plot(
                wavefunc.basis_labels / np.pi,
                ### adding in scaling factor for viewing purpose
                wavefunc.amplitudes*4.0 + energy_scaled, 
                '-',
                linewidth = 2.0
                )

            # plotting lines where wavefunctions asumptote
            self.subplots[plot_num].plot(
                wavefunc.basis_labels / np.pi,
                wavefunc.amplitudes*0.0 + energy_scaled, 
                'k-',
                linewidth = 0.5
                )

        # plot potential
        phi_range = np.linspace(-phi_limit*np.pi, phi_limit*np.pi, 500)
        potential = np.zeros(len(phi_range))
        for i in range(0, len(phi_range)):
            potential[i] = self.fluxonium.potential(phi_range[i])

        self.subplots[plot_num].plot(
            phi_range / np.pi,
            potential,
            'k-',
            linewidth = 2.0
            )
            
        self.subplots[plot_num].set_xlabel(r'$\varphi / \pi$')
        self.subplots[plot_num].set_ylabel('Energy (GHz)')
        self.subplots[plot_num].set_title('wave functions')


    ## method to plot Raman values against flux
    #
    def Raman_plot(self, plot_num = 0):
        # if method is called set to true
        self.Raman_bool = True
        # keep track of plot number
        self.Raman_plot_num = plot_num

        self.subplots[plot_num].cla()
        flux_points = 100 # number of flux points
        phi_ext = np.linspace(0.0, 1.0, flux_points)

        raman_range = np.zeros(flux_points)
        
        for i in range(0,flux_points):
            self.fluxonium.flux = phi_ext[i]
            raman_range[i] = np.abs(self.RamanVal(
                drive_freq = 0.05,
                beta_phi = 0.2,
                delta_lil = 0.005
                )) * 1e3 # MHz units

        self.subplots[plot_num].plot(phi_ext, raman_range, '-')

        self.subplots[plot_num].set_xlabel(r'$\varphi_{ext}/2\pi$')
        self.subplots[plot_num].set_ylabel('Raman Frequency (MHz)')
        self.subplots[plot_num].set_title(
            'Raman Frequency Plot')
        self.subplots[plot_num].set_yscale('log')


    ## method to calculte Raman frequency 
    # @param delta_lil: small detuning
    # @param beta_phi: beta parameter for coupling charactorization
    # 
    def RamanVal(self, 
            drive_freq = 0.05,
            delta_lil = 0.0,
            beta_phi = 0.2,
            ):
        # find dipole matrix
        g_mat = self.fluxonium.matrixelement_table(
            'n_operator',
            evals_count = self.Num_sum
            )

        g_mat = np.abs(g_mat)

        Energies = self.fluxonium.eigenvals(self.Num_sum)
        Energies = Energies - Energies[0]
        deltas = Energies - self.cav_freq

        # g_const is the coupling constant factor
        g_const = self.cav_freq * beta_phi * (2.0/np.pi)
        # drive is the drive factor
        drive = drive_freq / g_mat[0, 2]

        # arrays containing all transition values
        g_vals = g_const * g_mat[:, 1]
        omega_vals = drive * g_mat[0, :] 

        raman_vals = np.zeros(self.Num_sum)
        for i in range(1, self.Num_sum):
            raman_vals[i] = (
                g_vals[i] * omega_vals[i] / ((deltas[i] - delta_lil)*2) 
                )

        omega_raman = np.sum(raman_vals)

        return np.abs(omega_raman)

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
       self.fluxonium.EL = self.EL
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
            self.Flux_sweep_plot(plot_num = self.flux_sweep_num)
        if self.Raman_bool == True:
            self.Raman_plot(plot_num = self.Raman_plot_num)
        if self.g_mat_bool == True:
            self.g_mat_plot(plot_num = self.g_mat_num)
        if self.wavefunction_bool == True:
            self.wavefunction_plot(plot_num = self.wavefunction_num)

    ## method to show plots
    #
    def show(self):
        plt.show()


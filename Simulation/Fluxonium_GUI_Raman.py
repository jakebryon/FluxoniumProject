#!/usr/bin/python
# -*- coding: utf-8 -*-

#############################################################################
#############################################################################
#############################################################################
# Importing modules

"""
Created on Sun Feb 17 14:03:29 2019

@author: jacobbryon

Welcome to the code!

This code computes the energy level structure and the raman freqencies, amoung
other things, for fluxonium. It is a GUI type plot, so when it is run, pretty
much everything you want can be found just by pointing and clicking. It is 
important to note that this code runs a bit slow since there is alot going
on. 
Below I will layout a general idea of the things plotted in each figure, how
to 
use the intereface, and a few of the basic changes that you can make.

Figure 10: Fluxonium Energy Levels
This is the main figure and it sets the parameters that get plugged into the
other
figures. You will see two plots with some slidders underneath them. The
sliders
define, as they are labeled ont the left with the range of the bar and the 
paramter that it changes: the cavity frequency, EC, EJ, and EL. All the
energies 
given here are all in GHz. To change the paramter, just click on the bar where
you want to change it and it will update the plot (likely it will take several
moments to run). The left plot shows the transition energies of the fluxonium
levels, the 0th state energy has been set to zero. The brown horizontal line
shows the cavity frequency and the blue vertical line shows the flux value
that
is set for the other plots. To change the flux, the blue line, just click on
the
plot and it will move it to where you clicked, note that this takes longer
than
the other parameters do when changed. Occasionally I have had the plot freeze
up
on me, so you may need to quit the kernel and try again. The left plot just
shows
the potential and the plotted wave functions. As well the raman and g4 term
are
given for the particular inputs

Figure 20: Cavity Variation
Using the chosen paramters from the first figure, the left plot shows the
potential
with just the energy levels plotted over it. On the right is the absolute
value 
of the raman against varying cavity frequency. The blue, yellow, red, purple, 
vertical lines are each of the energy level transtions.

Figure 30: Couplings and Raman
Again, give the inputs from figure 10, we have 3 plots here all plotted
against 
change external flux. Top left are the 0th state couplings with the other
levels, 
similiarly the top right shows the same thing but with the 1st state. The
bottom
plot shows the raman and g4 terms as a function of external flux

Changing the Raman drive:
Around line 380, there are the initial values that are used in the plots.
Here 
under 'Raman Terms' are listed the drive_freq, which the the frequency that
the 
drive between the desired drive state and the excited state will be set too.
levelstart
gives which level will be driven, and delta_lil is the small detuning. So if
you
wanted to drive the zero-excited state transtion at 100MHz, choose the
following
drive_freq = .100
levelstart = 0 

I have gone through the code and tried to work out any bugs I have found and
cleaned it up, 
however I know that there are of course mistakes in here. If you have any
questions
or concerns, feel free to email me: jbryon@princeton.edu

"""

import numpy as np
import scipy.sparse as ssp
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.widgets import CheckButtons
import scipy.io
import matplotlib.gridspec as gridspec
from Fluxonium_BasicFunctions import *



###############################################################################
###############################################################################
###############################################################################
##########

plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 4

##########

# General constants

q_e = 1.602e-19  # electron charge in C
h = 6.626e-34  # Planck constant in Js
hbar = h / (2 * np.pi)  # Planck constant / 2p in Js
phi_0 = hbar / (2 * q_e)  # flux quantum / 2p in Wb

# beta_phi  = 0.2 # original given value for Beta

#beta_phi = 0.20  # hopeful value
beta_phi = 0.31  # FDH02 Simulation value Coupling Resonator
#beta_phi = 0.0446  # FDH02 Simulation value Readout Resonators
f_res = 7.8

g_phi = 0.125 * beta_phi * f_res #g_PHI = g_factor*charge_matrix element in GHz

tol = 1e-8  # diagonalization tolerance
keig = 8  # number of required ekets
N_grid = 251  # number of points in each direction (odd)
max_grid = 4 * np.pi  # maximum range of each coordinate
grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                       dtype='float64')
grid_x = keig
grid_y = 1
legend = []
scaling = 15


################################################################################
################################################################################
############# Creating functions for plotting and making GUI
###################################################################
######## setting up figures
# fig is the main plot that shows the energy levels and wave functions
# fig2 is just the energy levels and Raman vs cavity frequency
# fig_gs shows couplings and Raman and g4 as a function of phi_ext

axis_color = 'lightgoldenrodyellow'

### figures 1 and 2

fig = plt.figure(10, figsize=[12.0, 8.0])
fig2 = plt.figure(20, figsize=[12.0, 8.0])
plt.clf()
fig.clear()
fig2.clear()
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax3 = fig2.add_subplot(121)
ax4 = fig2.add_subplot(122)
legend = ['g', 'e1', 'e2', 'f']

#### figure 3, fig_gs

gs = gridspec.GridSpec(2, 2)
fig_gs = plt.figure(30, figsize=[12.0, 8.0])
fig_gs.clear()
ax_gs = fig_gs.add_subplot(gs[0, 0])
ax_gs2 = fig_gs.add_subplot(gs[0, 1])
ax_ramflux = fig_gs.add_subplot(gs[1, :])

#### figure 3, fig_gs

fig_Chi = plt.figure(40, figsize=[12.0, 8.0])
fig_Chi.clear()
ax_Chi1 = fig_Chi.add_subplot(211)
ax_Chi2 = fig_Chi.add_subplot(212)


# Adjust the subplots region to leave some space for the sliders and buttons

fig.subplots_adjust(left=0.1, bottom=0.35)
fig2.subplots_adjust(left=0.1, bottom=0.35)

##############################################################################
######################################### intial values for plots
#### START

phi_range = [0, 1]
N_phi = 51
t = np.linspace(phi_range[0] * 2 * np.pi, phi_range[1] * 2 * np.pi,
                N_phi) / 2 / np.pi

# fluxonium and cavity
### All energies are in GHz

#### testing Shuster long T1 parameters
#EC = 0.24
#EJ = 8.11
#EL = 0.46
#cavfreq = 4.95

### Testing parameters from Vlad Thesis
#EC = 2.500
#EJ = 8.900
#EL = .525
#cavfreq = 8.175

### Testing parameters from Anjali Qubit
#EC = 0.44
#EJ = 3.85
#EL = 0.25
#cavfreq = 7.50


#### FDH01 Qubit 1 parameters
#EC = 1.15
#EJ = 11.25
#EL = 0.96
#cavfreq = 8.74

#### FDH01 Qubit 2 parameters
#EC = 1.19
#EJ = 10.78
#EL = 0.96
##cavfreq = 7.58
##cavfreq = 7.0 ### FDH02 Replica Read Out Resonator 1
#cavfreq = 8.0 ### FDH02 Replica Read Out Resonator 1

#########################################
##  new parameters for FDH02  2
EC = 1.16
EJ = 8.34
EL = 0.45
cavfreq = 6.0
#cavfreq = 5.5 ### reaadout 1
#cavfreq = 6.5 ### readout 2
#######################################



drive_freq = .05
delta_lil = .05
levelstart = 0  # level used for drive

# starting display phi ext value

xval = 0.6
phi_ext = xval * 2 * np.pi

### ploting lines for showing cuts on plots

vert = ax.axvline(x=xval)

# line for the varying phi ext plot

cavline = ax.axhline(y=cavfreq, color='olive')

# line for wave function plot

cavhorz = ax2.axhline(y=cavfreq, color='olive')
ground1 = ax4.axvline(x=0)
excit1 = ax4.axvline(x=0)
excit2 = ax4.axvline(x=0)


####################################################################
###########################################################################

def fluxoniumLevels(
    EC,
    EJ,
    EL,
    cavfreq,
    ):
    global vert
    global cavhorz
    global ground1
    global excit1
    global excit2

    #### drawing fluxonium potential and energies, 
    #### with energies against ext flux
    #### drawing ax, energy levels against flux

    ax.clear()
    phi_ext_s = np.linspace(phi_range[0] * 2 * np.pi, phi_range[1] * 2
                            * np.pi, N_phi)
    (evals_s, ekets_s) = flux_sweep(phi_ext_s, EC, EJ, EL)
    [line1, ] = ax.plot(
        t,
        (evals_s.T - evals_s[:, 0])[0],
        linewidth=2,
        alpha=0.5,
        color='blue',
        marker='o',
        )
    [line2, ] = ax.plot(
        t,
        (evals_s.T - evals_s[:, 0])[1],
        linewidth=2,
        alpha=0.5,
        color='orange',
        marker='o',
        )
    [line3, ] = ax.plot(
        t,
        (evals_s.T - evals_s[:, 0])[2],
        linewidth=2,
        alpha=0.5,
        color='green',
        marker='o',
        )
    [line4, ] = ax.plot(
        t,
        (evals_s.T - evals_s[:, 0])[3],
        linewidth=2,
        alpha=0.5,
        color='red',
        marker='o',
        )
    [line5, ] = ax.plot(
        t,
        (evals_s.T - evals_s[:, 0])[4],
        linewidth=2,
        alpha=0.5,
        color='purple',
        marker='o',
        )
    vert = ax.axvline(x=xval)

    # ### final plot edits

    legend = ['g1', 'g2', 'e1', 'e2', 'e3']

    cavline = ax.axhline(y=cavfreq, color='olive')

    ax.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax.set_ylabel('E - E(0) (GHz)')
    ax.legend(
        legend,
        bbox_to_anchor=(-0.25, 0., 0.20, -0.75),
        loc=3,
        ncol=2,
        mode='expand',
        borderaxespad=0.,
        )
    ax.set_title('Energy Levels vs. ' + r'$\varphi_{ext}/2\pi$')

    # ###### drawing ax2

    ax2.clear()
    phi_ext = xval * 2 * np.pi

    keig = 10
    keig2 = keig  # number of required ekets
    N_grid = 251  # number of points in each direction (odd)
    max_grid = 4 * np.pi  # maximum range of each coordinate

    H_op = H_inductive(phi_ext, EC, EJ, EL)
    U = H_potential(phi_ext, EJ, EL)
    (evals_s, ekets_s) = diagonalize_new(H_op, keig)

    grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                           dtype='float64')

    grid_x = 5
    grid_y = 1
    legend = []
    scaling = 15

    minlist = np.zeros(grid_x)
    maxlist = np.zeros(grid_x)
    evals_s = evals_s - evals_s[0]

    for i in range(grid_x):
        for j in range(grid_y):
            idx = grid_y * i + j
            ax2.plot(grid_pts / np.pi, scaling * ekets_s[:, idx]
                     + evals_s[idx], alpha=0.8)

            # ## finding y limits

            minlist[i] = np.min(scaling * ekets_s[:, idx]
                                + evals_s[idx])
            maxlist[i] = np.max(scaling * ekets_s[:, idx]
                                + evals_s[idx])

            # ##

            legend.append('Eigval = %.3f' % (evals_s[idx] - evals_s[0]))
            ax2.legend(legend, bbox_to_anchor=(0.75, -.12), loc=2,
                       borderaxespad=0.)
            ax2.set_xlabel(r'$\varphi/\pi$')

    ylow = np.min(minlist) - 10
    yhigh = np.max(maxlist) + 10

    cavhorz = ax2.axhline(y=cavfreq, color='olive')
    legend.append('Cavity Frequency')
    ax2.legend(legend, bbox_to_anchor=(0.75, -.12), loc=2,
               borderaxespad=0.)

    ax2.set_ylim([ylow, yhigh])
    ax2.plot(grid_pts / np.pi, U, '-k')
    ax2.set_ylabel(r'$\omega/2\pi$ [GHz]')
    ax2.set_title(r'Wavefunctions for $\varphi_{ext}/2\pi$ = %.4f'
                  % (phi_ext / 2 / np.pi))

    OmgRaman = RamanVal(
        EC,
        EJ,
        EL,
        cavfreq,
        phi_ext,
        levelstart,
        )

    gfourth = gfourthVal(EC, EJ, EL, cavfreq, phi_ext)

    ChiVal = ChiShift(
            EC,
            EJ,
            EL,
            cavfreq,
            phi_ext,
            level = levelstart,
            )

    ax2.text(2.2, -7.5, r'$\Omega_{Raman}$ = %.4f' % (OmgRaman * 1e3)
             + ' MHz')
    ax2.text(2.2, -9.0, r'$\chi$ = %.4f' % (ChiVal *1e3)
            + ' MHz')
    ax2.text(2.2, -10.5, r'$g^{4}$ term = %.4f' % (gfourth * 1e3)
             + ' MHz')

    fig.suptitle(r'Fluxonium Energy Levels, Raman Level = %i'
                 % levelstart)


##########################################################################
###########################################################################

def fluxoniumRamanCavity(
    EC,
    EJ,
    EL,
    cavfreq,
    ):
    global vert
    global cavhorz
    global ground1
    global excit1
    global excit2

    ######### plotting  Raman against Cavity Frequency

    phi_ext = xval * 2 * np.pi

    keig = 5
    grid_x = keig
    grid_y = 1
    legend = []
    scaling = 15
    keig2 = keig

    H_op = H_inductive(phi_ext, EC, EJ, EL)
    U = H_potential(phi_ext, EJ, EL)
    (evals_s, ekets_s) = diagonalize_new(H_op, keig2)
    evals_s = evals_s - evals_s[0]  # setting zero energy

    minlist = np.zeros(grid_x)
    maxlist = np.zeros(grid_x)

    ##### finding the dipole matrix

    (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)
    gconst = cavfreq * beta_phi * (2 / np.pi)
    g_matrix = g_matrix * gconst

    ##########################################################################

    grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                           dtype='float64')

    grid_x = keig
    grid_y = 1
    legend = []
    scaling = 15

    ax3.clear()

    # ######### creating another plot for only the level structure

    eng_levels = np.ones((len(grid_pts), len(evals_s)))
    i = 0
    for i in range(grid_x):
        for j in range(grid_y):
            idx = i
            eng_levels[:, idx] = evals_s[idx] * eng_levels[:, idx]
            ax3.plot(grid_pts / np.pi, eng_levels[:, idx], alpha=0.8)

            # ## finding y limits

            minlist[i] = np.min(scaling * ekets_s[:, idx]
                                + evals_s[idx])
            maxlist[i] = np.max(scaling * ekets_s[:, idx]
                                + evals_s[idx])

            # ##

            legend.append('Eigval = %.3f' % (evals_s[idx] - evals_s[0]))
            ax3.legend(legend)
            ax3.set_xlabel(r'$\varphi/\pi$')

    ylow = np.min(minlist) - 10
    yhigh = np.max(maxlist) + 10

    ax3.legend(legend, bbox_to_anchor=(0., -.12), loc=2,
               borderaxespad=0.)

    ax3.set_ylim([ylow, yhigh])
    ax3.plot(grid_pts / np.pi, U, '-k')
    ax3.set_ylabel(r'$\omega/2\pi$ [GHz]')
    ax3.set_title(r'EC = %.2f' % EC + r'; EC = %.2f' % EJ
                  + r'; EC = %.2f' % EL)

    # ##### finding the raman transition and plotting

    ax4.clear()
    legend = []

    # #### finding the dipole matrix

    (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)

    numomg = 300
    CavList = np.linspace(0, 20, numomg)
    RamanList = np.ones(numomg)

    for i in range(0, numomg):
        RamanList[i] = RamanVal(
            EC,
            EJ,
            EL,
            CavList[i],
            phi_ext,
            numstart=levelstart,
            )

    ground0 = ax4.axvline(x=evals_s[0], color='blue')
    ground1 = ax4.axvline(x=evals_s[1], color='orange')
    excit1 = ax4.axvline(x=evals_s[2], color='green')
    excit2 = ax4.axvline(x=evals_s[3], color='red')
    excit3 = ax4.axvline(x=evals_s[4], color='purple')

    ax4.plot(CavList, np.abs(RamanList) * 1e3, color='olive')
    ax4.legend([
        'g1',
        'g2',
        'e1',
        'e2',
        'e3',
        'Raman',
        ], bbox_to_anchor=(0., -.12), loc=2, borderaxespad=0.)
    ax4.set_xlabel('Cavity Frequency (GHz)')
    ax4.set_ylabel('$\Omega_{Raman}$ (MHz)')
    ax4.set_title('Abs Raman Frequency')
    ax4.set_yscale('log')

    fig2.suptitle(r'Cavity Variation $\varphi_{ext}/2\pi$ = %.4f'
                  % (phi_ext / 2 / np.pi) + r', Raman Level = %i'
                  % levelstart)


    # ##########################################################################

###########################################################################

def fluxoniumcouplingPlots(
    EC,
    EJ,
    EL,
    cavfreq,
    ):
    global vert
    global cavhorz
    global ground1
    global excit1
    global excit2

    # #### ploting matrix elements

    keig2 = 10

    # ### for loop parameters

    range_limit = 200
    gmat_list = np.zeros((keig2, keig2, range_limit))
    phiext_list = np.linspace(0, 2 * np.pi, range_limit)
    RamanList = np.ones(range_limit)
    ChiList = np.ones(range_limit)


    # ## finding the raman transition

    for ind in range(0, range_limit):
        phi_ext = phiext_list[ind]
        H_op = H_inductive(phi_ext, EC, EJ, EL)
        (evals_s, ekets_s) = diagonalize_new(H_op, keig2)
        evals_s = evals_s - evals_s[0]  # setting zero energy

        # finding the dipole matrix

        gconst = cavfreq * beta_phi * (2 / np.pi)
        (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)
        gmat_list[:, :, ind] = gconst * g_matrix

        RamanList[ind] = RamanVal(
            EC,
            EJ,
            EL,
            cavfreq,
            phi_ext,
            numstart=levelstart,
            )

    # #### finding g4 values

    g4totalList = np.zeros(range_limit)

    for ind in range(0, range_limit):
        phi_ext = phiext_list[ind]
        g4totalList[ind] = gfourthVal(EC, EJ, EL, cavfreq, phi_ext)

    # ########

    # ### couplings

    g0_list = gmat_list[:, 0, :] * 1e3
    g1_list = gmat_list[:, 1, :] * 1e3

    # ## plotting

    leg_list = np.arange(keig2)
    legend = leg_list

    ax_gs.clear()
    ax_gs2.clear()
    ax_ramflux.clear()
    phiext_list = phiext_list / (2 * np.pi)

    ax_gs.plot(phiext_list, np.abs(g0_list[0]))
    ax_gs.plot(phiext_list, np.abs(g0_list[1]))
    ax_gs.plot(phiext_list, np.abs(g0_list[2]))
    ax_gs.plot(phiext_list, np.abs(g0_list[3]))
    ax_gs.plot(phiext_list, np.abs(g0_list[4]))
    ax_gs.plot(phiext_list, np.abs(g0_list[5]))
    ax_gs.plot(phiext_list, np.abs(g0_list[6]))
    ax_gs.plot(phiext_list, np.abs(g0_list[7]))
    ax_gs.plot(phiext_list, np.abs(g0_list[8]))
    ax_gs.plot(phiext_list, np.abs(g0_list[9]))
    ax_gs.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax_gs.set_ylabel('absolute value g value MHz')
    ax_gs.legend(legend, bbox_to_anchor=(-0.3, 1.0), loc=2,
                 borderaxespad=0.)
    ax_gs.set_title('g0 Matrix Elements')

    ax_gs2.plot(phiext_list, np.abs(g1_list[0]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[1]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[2]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[3]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[4]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[5]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[6]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[7]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[8]))
    ax_gs2.plot(phiext_list, np.abs(g1_list[9]))
    ax_gs2.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax_gs2.set_ylabel('absolute value g value MHz')
    ax_gs2.legend(legend, bbox_to_anchor=(1.1, 1.0), loc=2,
                  borderaxespad=0.)
    ax_gs2.set_title('g1 Matrix Elements')

    # ## plotting

    ax_ramflux.plot(phiext_list, np.abs(RamanList) * 1e3)
    ax_ramflux.plot(phiext_list, np.abs(g4totalList) * 1e3)
    #ax_ramflux.plot(phiext_list, np.abs(ChiList) * 1e3)
    ax_ramflux.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax_ramflux.set_ylabel('absolute value $\Omega_{Raman}$ (MHz)')
    ax_ramflux.set_title(r'Raman for Cavity %.4f' % cavfreq + 'GHz')
    ax_ramflux.legend(['Raman', 'g4', 'Chi'])
    ax_ramflux.set_yscale('log')

    fig_gs.suptitle(r' Couplings and Raman, EC = %.2f' % EC
                    + r'; EJ = %.2f' % EJ + r'; EL = %.2f' % EL
                    + r', Raman Level = %i' % levelstart)


###########################################################################
def fluxoniumChiShift(
    EC,
    EJ,
    EL,
    cavfreq,
    ):
    global vert
    global cavhorz
    global ground1
    global excit1
    global excit2

    keig2 = 10
    ChiNum = 5 #number of levels to look at
    #### for loop parameters

    range_limit = 200
    phiext_list = np.linspace(0, 2 * np.pi, range_limit)
    ChiList = np.ones((range_limit,ChiNum))

    # ## finding the raman transition

    for ind in range(0, range_limit):
        phi_ext = phiext_list[ind]
        H_op = H_inductive(phi_ext, EC, EJ, EL)
        (evals_s, ekets_s) = diagonalize_new(H_op, keig2)
        evals_s = evals_s - evals_s[0]  # setting zero energy

        # finding the dipole matrix

        gconst = cavfreq * beta_phi * (2 / np.pi)
        (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)

        for ind2 in range(0, ChiNum):
            ChiList[ind, ind2] = ChiShift(
                EC,
                EJ,
                EL,
                cavfreq,
                phi_ext,
                level = ind2
                )

    ## plotting
    ax_Chi1.clear()
    ax_Chi2.clear()
    phiext_list = phiext_list / (2 * np.pi)
    ### testing abs value
    #ChiList = np.abs(ChiList)


    ##### inidividual level 
    legend1 = ['0', '1', '2', '3', '4'] 

    ax_Chi1.plot(phiext_list, ChiList[:, 0]*1e3)
    ax_Chi1.plot(phiext_list, ChiList[:, 1]*1e3)
    ax_Chi1.plot(phiext_list, ChiList[:, 2]*1e3)
    ax_Chi1.plot(phiext_list, ChiList[:, 3]*1e3)

    ax_Chi1.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax_Chi1.set_ylabel('Energy MHz')
    ax_Chi1.legend(legend1)
    ax_Chi1.set_title('Individual Level Shifts')
    #ax_Chi1.set_yscale('log')
    #ax_Chi1.set_ylim([-1000, 1000])

    ###### diff plot
    legend2 = ['0', '1', '2', '3', '4'] 
    
    ax_Chi2.plot(phiext_list, (ChiList[:, 0] - ChiList[:,0])*1e3)
    ax_Chi2.plot(phiext_list, (ChiList[:, 1] - ChiList[:,0])*1e3)
    ax_Chi2.plot(phiext_list, (ChiList[:, 2] - ChiList[:,0])*1e3)
    ax_Chi2.plot(phiext_list, (ChiList[:, 3] - ChiList[:,0])*1e3)

    ax_Chi2.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax_Chi2.set_ylabel('Energy MHz')
    ax_Chi2.legend(legend2)
    ax_Chi2.set_title('Level Shifts from 0')

    ax_ramflux.set_yscale('log')

    fig_Chi.suptitle(r' Chi Shifts, cav = %.2f ' %cavfreq + r'; EC = %.2f' % EC
                    + r'; EJ = %.2f' % EJ + r'; EL = %.2f' % EL
                    )


###########################################################################


##################################
### plotting function

def plot(
    EC,
    EJ,
    EL,
    cavfreq,
    ramlevel=0,
    ):
    global vert
    global cavhorz
    global ground1
    global excit1
    global excit2

    # ########################################################################

    fluxoniumLevels(EC, EJ, EL, cavfreq)
    # commenting out secondary plots
    #fluxoniumRamanCavity(EC, EJ, EL, cavfreq)
    fluxoniumcouplingPlots(EC, EJ, EL, cavfreq)
    fluxoniumChiShift(EC, EJ, EL, cavfreq)


   #########################################################################

##################################################################
###################################################################

# Adding sliders for tweaking the parameters

Cav_slider_ax = fig.add_axes([0.1, 0.20, 0.65, 0.03],
                             facecolor=axis_color)
Cav_slider = Slider(Cav_slider_ax, r'$\omega_{cav}$ [3:10]', 3.0, 10.0,
                    valinit=cavfreq)
#######
EC_slider_ax = fig.add_axes([0.1, 0.15, 0.65, 0.03],
                            facecolor=axis_color)
EC_slider = Slider(EC_slider_ax, 'EC [0.2:3]', 0.20, 3.0, valinit=EC)
#######
EJ_slider_ax = fig.add_axes([0.1, 0.1, 0.65, 0.03],
                            facecolor=axis_color)
EJ_slider = Slider(EJ_slider_ax, 'EJ [3:12]', 3.0, 12.0, valinit=EJ)
#######
EL_slider_ax = fig.add_axes([0.1, .05, 0.65, 0.03],
                            facecolor=axis_color)
EL_slider = Slider(EL_slider_ax, 'EL [0.1:1]', 0.1, 1.0, valinit=EL)


###################################################################
#### creating slider change action

def sliders_on_changed(val):
    global EC
    global EJ
    global EL
    global cavfreq
    plot(EC_slider.val, EJ_slider.val, EL_slider.val, Cav_slider.val)
    EC = EC_slider.val
    EJ = EJ_slider.val
    EL = EL_slider.val
    cavfreq = Cav_slider.val

    fig.canvas.draw_idle()
    fig2.canvas.draw_idle()
    fig_gs.canvas.draw_idle()
    fig_Chi.canvas.draw_idle()

###################################################################
###################################################################
# Add a button for resetting the parameters

reset_button_ax = fig.add_axes([0.85, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color,
                      hovercolor='0.975')


def reset_button_on_clicked(mouse_event):
    global xval
    global vert
    global phi_ext
    phi_ext = xval * 2 * np.pi
    EJ_slider.reset()
    EC_slider.reset()
    EL_slider.reset()
    Cav_slider.reset()

    # ##########################3 resetting ax2

    xval = 0.6
    vert.set_xdata(xval)
    plot(EC, EJ, EL, cavfreq)


###################################################################
###### adding clicking mechanism

def onclick(event):
    global xval
    global vert
    global cavline
    global phi_ext

    if event.inaxes in [ax]:
        xval = event.xdata
        phi_ext = 2 * np.pi * xval
        vert.set_xdata(xval)
        cavline.set_ydata(cavfreq)
        plot(EC, EJ, EL, cavfreq)
        fig.canvas.draw_idle()
        fig2.canvas.draw_idle()
        fig_gs.canvas.draw_idle()
        fig_Chi.canvas.draw_idle()

    else:
        return


###################################################################
###  note that the original figure creation is above, to turn off
### figure entirely, must comment out figure
###### turning on sliders

Cav_slider.on_changed(sliders_on_changed)
EC_slider.on_changed(sliders_on_changed)
EJ_slider.on_changed(sliders_on_changed)
EL_slider.on_changed(sliders_on_changed)

############ turning on reset button

reset_button.on_clicked(reset_button_on_clicked)

### turn on clicking

cid = fig.canvas.mpl_connect('button_press_event', onclick)

### intial plot

plot(EC, EJ, EL, cavfreq)
plt.show()

#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
################################################################################
################################################################################
# Importing modules

"""
Created on Sun Feb 17 14:03:29 2019

@author: jacobbryon
"""

import time
import math
import numpy as np
import scipy.sparse as ssp
import scipy.sparse.linalg
from scipy.sparse import diags
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from IPython.display import display, Image, Markdown
from scipy.io import loadmat
import timeit
import scipy.io
from mpl_toolkits.mplot3d import axes3d
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

beta_phi = 0.2
f_res = 7.8

g_phi = 0.125 * beta_phi * f_res  # g_PHI = g_factor * charge_matrix element in GHz

tol = 1e-8  # diagonalization tolerance
keig = 5  # number of required ekets
N_grid = 251  # number of points in each direction (odd)
max_grid = 4 * np.pi  # maximum range of each coordinate
grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                       dtype='float64')
grid_x = keig
grid_y = 1
legend = []
scaling = 15


###############################################################################
### Creating functions and making GUI
###################################################################
######### START
######## setting up figure

axis_color = 'lightgoldenrodyellow'

fig = plt.figure(10, figsize=[12.0, 8.0])
plt.clf()
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
legend = ['g', 'e1', 'e2', 'f']

# Adjust the subplots region to leave some space for the sliders and buttons

fig.subplots_adjust(left=0.1, bottom=0.35)

##### intial values for plots

phi_range = [0, 1]
N_phi = 51
t = np.linspace(phi_range[0] * 2 * np.pi, phi_range[1] * 2 * np.pi,
                N_phi) / 2 / np.pi

# fluxonium and cavity
### All energies are in GHz

#### testing Shuster long T1 parameters
#EC = 0.24
#EJ = 8.11
##EL = 0.46
#EL = 0.0
#cavfreq = 4.95

### Testing parameters from Vlad Thesis
#EC = 2.500
#EJ = 8.900
#EL = .525
#cavfreq = 8.175

#### FDH01 Qubit 1 parameters
#EC = 1.15
#EJ = 11.25
#EL = 0.96
#cavfreq = 8.74

### FDH01 Qubit 2 parameters
EC = 1.19
EJ = 10.78
EL = 0.96
cavfreq = 7.58

xval = 0.342
vert = ax.axvline(x=xval)
### adding cavity for fining ChiShift
levelstart = 0

##################################################################
##################################
### plotting function

def plot(EC, EJ, EL):
    global vert

    # ################
    # #### drawing ax

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

#    vert.set_xdata(x=xval)

    vert = ax.axvline(x=xval)

    # ### final plot edits

    legend = ['g1', 'g2', 'e1', 'e2', 'e3']

    ax.set_xlabel(r'$\varphi_{ext}/2\pi$')
    ax.set_ylabel('E - E(0) (GHz)')
    ax.legend(
        legend,
        bbox_to_anchor=(0., -0.28, 1., -0.5),
        loc=3,
        ncol=2,
        mode='expand',
        borderaxespad=0.,
        )
    ax.set_title('Energy Levels vs. ' + r'$\varphi_{ext}/2\pi$')

    # ###### drawing ax2

    ax2.clear()
    phi_ext = xval * 2 * np.pi

    # ######## All energies are in GHz

    keig = 5  # number of required ekets
    N_grid = 251  # number of points in each direction (odd)
    max_grid = 4 * np.pi  # maximum range of each coordinate

    H_op = H_inductive(phi_ext, EC, EJ, EL)
    U = H_potential(phi_ext, EJ, EL)
    (evals_s, ekets_s) = diagonalize_new(H_op, keig)

    grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                           dtype='float64')

    #    plt.figure(num = 2, figsize=(6,6))

    grid_x = keig
    grid_y = 1
    legend = []
    scaling = 15

    minlist = np.zeros(grid_x)
    maxlist = np.zeros(grid_x)
    evals_s = evals_s - evals_s[0]

    # horz1 = ax2.axhline()

    for i in range(grid_x):
        for j in range(grid_y):
            idx = grid_y * i + j
            ax2.plot(grid_pts / np.pi, scaling * ekets_s[:, idx]
                     + evals_s[idx], alpha=0.8)

            # ####printing horizontal to show energy of level
    #        horz1.set_ydata(evals_s[idx])
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

    ax2.set_ylim([ylow, yhigh])
    ax2.plot(grid_pts / np.pi, U, '-k')
    ax2.set_ylabel(r'$\omega/2\pi$ [GHz]')
    ax2.set_title(r'Wavefunctions for $\varphi_{ext}/2\pi$ = %.4f'
                  % (phi_ext / 2 / np.pi))

###################################################################
###################################################################
# Adding sliders for tweaking the parameters

Cav_slider_ax = fig.add_axes([0.1, 0.20, 0.65, 0.03],
                             facecolor=axis_color)
Cav_slider = Slider(Cav_slider_ax, r'$\omega_{cav}$ [5:10]', 5.0, 10.0,
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
    plot(EC_slider.val, EJ_slider.val, EL_slider.val)

    fig.canvas.draw_idle()


###################################################################

###################################################################
# Add a button for resetting the parameters

reset_button_ax = fig.add_axes([0.85, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color,
                      hovercolor='0.975')


def reset_button_on_clicked(mouse_event):
    global xval
    global vert
    EJ_slider.reset()
    EC_slider.reset()
    EL_slider.reset()

    # ##########################3 resetting ax2

    xval = 0.75
    vert.set_xdata(xval)
    plot(EC, EJ, EL)


###################################################################
###### adding clicking mechanism

def onclick(event):
    global xval
    global vert
    global phi_ext
    global EC
    global EJ
    global EL
 
    if event.inaxes in [ax]:
        xval = event.xdata
        phi_ext = 2*np.pi * xval
        vert.set_xdata(xval)
        plot(EC, EJ, EL)
        fig.canvas.draw_idle()
    else:
        return


###################################################################
###  note that the original figure creation is above, to turn off
### figure entirely, must comment out figure
###### turning on sliders

EC_slider.on_changed(sliders_on_changed)
EJ_slider.on_changed(sliders_on_changed)
EL_slider.on_changed(sliders_on_changed)

############ turning on reset button

reset_button.on_clicked(reset_button_on_clicked)

### turn on clicking

cid = fig.canvas.mpl_connect('button_press_event', onclick)

### intial plot

plot(EC, EJ, EL)
plt.suptitle('Fluxonium Energy Levels')
plt.show()

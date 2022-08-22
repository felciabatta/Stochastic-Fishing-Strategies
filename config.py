#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 18:31:00 2022

@author: felixdubicki-piper
"""

import math as m
import matplotlib as mpl

# ================= ==========================================================
#  RESULTS OPTIONS
# ================= ==========================================================

N = 1  # No. of iterations (to produce N results, higher -> more accurate)
saveResults = False
manyBoatPaths = False  # new path each iteration (False -> faster results)


# =========== =========================== ====================================
#  VARIABLES   [stuff you can mess with]
# =========== =========================== ====================================

# NOTE: (~1000**2 m**2 mean lake size, but slow to run, min 300**2)
A = 400**2     # lake area, m**2
L = m.sqrt(A)  # lake length, m

# NOTE: (~4 m/s typical fishing boat speed)
b = 4     # boat speed, m/s
bRad = 5  # boat catching radius, m (5 m for fishing line extension)

# NOTE: (0.35 is 'critical' speed for average Perch, Burst speed = 1.5m/s)
f = 0.35  # fish speed, m/s
fRad = 1  # fish shoal radius, m (see research ~1 default)

# ~~~~~~~~~~~~~~~~~~~~~
#  STARTING CONDITIONS  {angle and coordinate}
# ~~~~~~~~~~~~~~~~~~~~~

# boat
bL = L              # restrict boat to square of length bL, m (default = L)
bCoord0 = "random"  # [(L*0.99),L/2]  # "random" or [x, y]
bAngle0 = "random"  # "random" or float value

# fish
fCoord0 = "random"  # "random" or [x, y]
fAngle0 = "random"  # "random" or float value


# =========== =========================== ====================================
#  CONSTANTS   [try not to change these]
# =========== =========================== ====================================

maxT = 2*(fRad + bRad)/b - 0.1  # maxT*b <  2(fRad+bRad)
T = maxT  # time interval, set so    T*b <= 2(fRad+bRad), default = maxT

nrmDst = True  # normal dist vs flat dist for random direction
p = 0.75       # prob of changing direction each 'instant' (flat dist only)

# ~~~~~~~~~~~~~~~~~~~~~~~~
#  PROBABILITY PARAMATERS  {determines random change in angle}
# ~~~~~~~~~~~~~~~~~~~~~~~~
'''
The values here are the standard deviations for the normal distribution that
determines the change in angle every 'instant'.

Since the output of the norm dist. is mulplied by π, then the actual stn devs,
are π*bS and π*fS. i.e. bS and fS are a fraction of π.

So for S = 0.5 would mean an actual stn dev of π/2, so in a norm dist, that
means that 64% of the time, the angle change is between -π/2 and π/2.

For very high stn devs i.e. S > 1, the random change in angle becomes approx.
even distributed.

For very low stnd devs i.e. almost 0, motion is approx. linear.

Note that a smaller time period T results in more frequent changes in angle,
so to reproduce similar movement for all possible T, the stn dev should be
roughly proportional to T. i.e. for frequent changes in angle, we want the
possible change in angle to be smaller, so that it balances out and vice versa.
'''

# boat (0 for linear, 0.02+ for high std)
bS = 0.2  # std of boat angle

# fish (best to leave this as default)
proportionalS = 0.05*T+0.05  # std proportional to T
fS = proportionalS  # std of fish angle

# ================== =========================================================
#  PLOT & ANIMATION
# ================== =========================================================

plot = False  # produce plot of each iteration (will slow program down)
scatter = False  # scatter or line plot
if plot:
    title = input("Enter a title: ")
else:
    title = 'Title'
animated = False  # plot every instant separately (MIGHT CRASH PC)
ps = 10  # playback speed, 1 = realtime, max ~ 15 (MIGHT CRASH PC)

flabel = '$\sigma_{fish}$ = '+str(fS)
blabel = '$\sigma_{fish}$ = '+str(bS)

# ======= ====================================================================
#  FONTS
# ======= ====================================================================

mpl.rcdefaults()

mpl.rcParams['font.family'] = ['serif']

mpl.rcParams['font.serif'] = ['CMU Serif']+mpl.rcParamsDefault['font.serif']

mpl.rcParams['font.sans-serif'] = ['CMU Sans Serif'] + \
    mpl.rcParamsDefault['font.sans-serif']

mpl.rcParams['mathtext.fontset'] = 'cm'

# ======== ===================================================================
#  COLORS
# ======== ===================================================================

sweetpink = '#ff4ab9'
calmblue = '#00a8e3'

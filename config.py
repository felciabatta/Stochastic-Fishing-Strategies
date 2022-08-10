#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 18:31:00 2022

@author: felixdubicki-piper
"""
#CONFIG#

import math as m

# Results Options ---------------------------------------------------------------

N = 1 # No. of iterations (to produce N results, higher -> more accurate)
saveResults = False
manyBoatPaths = False # New boat path every iteration (False -> faster results)


# Variables (stuff you can mess with) -------------------------------------------


A = 100**2 # lake area, m^2 (~1000^2 m^2 average lake size, but takes long to run, min 300**2)
L = m.sqrt(A) # lake length, m

b = 4 # boat speed, m/s (~ 4 typical fishing boat speed)
bRad = 5 # boat catching radius, m (5m for fishing line extension)

f = 0.35 # fish speed, m/s (0.35 is 'critical' speed for average Perch, Burst speed = 1.5m/s)
fRad = 1 # fish shoal radius, m (see research ~ 1 default)



    # Starting conditions (angle and coordinate)

# boat
bL = L # restricts boat to smaller square, length m (default = L)
bCoord0 = "random" #[(L*0.99),L/2] # "random" or [x, y]
bAngle0 = "random" # "random" or float value

# fish
fCoord0 = "random" # "random" or [x, y]
fAngle0 = "random" # "random" or float value



# Constants (Try not to Change these) -------------------------------------------

maxT = 2*(fRad + bRad)/b - 0.1
T = maxT # time interval: set as maxT or such that T*b <= 2(fRad+bRad)

nrmDst = True # normal dist vs flat dist for random direction
p = 0.75 # probabilty of changing direction each 'instant' (Flat dist only)



    # Probabilty Parameters (determines random change in angle)

'''
The values here are the standard deviations for the normal distribution that
determines the change in angle every 'instant'.

Since the output of the norm dist. is mulplied by π, then the actual stn devs,
are π*bS and π*fS. i.e. bS and fS are a fraction of π.

So for S = 0.5 would mean an actual stn dev of π/2, so in a norm dist, that means
that 64% of the time, the angle change is between -π/2 and π/2.

For very high stn devs i.e. S > 1, the random change in angle becomes approx.
even distributed.

For very low stnd devs i.e. almost 0, motion is approx. linear.

Note that a smaller time period T results in more frequent changes in angle,
so to reproduce similar movement for all possible T, the stn dev should be
roughly proportional to T. i.e. for frequent changes in angle, we want the possible
change in angle to be smaller, so that it balances out and vice versa.
'''

# boat
bS = 0.2 # variation in boat direction (0 for linear, 0.02+ for high variation)

# fish (best to leave this as default)
proportionalS = 0.05*T+0.05
fS = proportionalS # variation in fish direction (proportional to T)



# Visuals - PLOT & ANIMATION ----------------------------------------------------

plot = True # produce plot of each iteration (will slow program down)
scatter = False # scatter or line plot
if plot:
    title = input("Enter a title: ")
else:
    title = None
animated = False # plot every instant separately (MIGHT CRASH PC)
ps = 10 # playback speed, 1 = realtime, max ~ 15 (MIGHT CRASH PC)

flabel = '$\sigma_{fish}$ = '+str(fS)
blabel = '$\sigma_{fish}$ = '+str(bS)

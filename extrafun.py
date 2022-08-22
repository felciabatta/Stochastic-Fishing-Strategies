#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 17:49:58 2022

@author: felixdubicki-piper
"""

import math as m
from simulation import *
from config import *


def plotFishTrajectory(tLim=18000, plotNum=1, σ=fS, speed=f):
    for i in range(plotNum):
        fishPath = createBoatPath(coord0=fCoord0,
                                  angle0=fAngle0,
                                  stnDev=σ,
                                  speed=speed)
        plotBoat(fishPath, tLim=tLim, col='b', label=f"Fish\nt={tLim}s",
                 title="Fish Trajectory, "+f"$\sigma={σ},$ $f={speed}$ $m/s$")


def plotBoatSpiral(tLim=1800, Q=0.6, period=450):
    boatPath = circularMotion(constRad=False,
                              RATIO=(1-Q),
                              radTimePeriod=period,
                              speed=b)
    plotBoat(boatPath, tLim=tLim, col='r', label=f"Boat\nt={tLim}s",
             title="Spiral, $T_{period}=$"+str(period)+"s, $Ratio=$"+str(Q))


def plotBoatTrajectory(tLim=1800, xy0=[0.0, 0.0], α0=m.pi/3, σ=0,
                       title=title, label=blabel):
    boatPath = createBoatPath(coord0=xy0,
                              angle0=α0,
                              stnDev=σ,
                              speed=b)
    plotBoat(boatPath, tLim=tLim, col='r', label=label, title=title)

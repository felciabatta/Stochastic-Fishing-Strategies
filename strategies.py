#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 17:28:44 2022

@author: felixdubicki-piper
"""
#STRATEGIES#

import math as m
import time
from simulation import *
from config import *

# run for different starting angles (betwen 45º and 90º exclusive)
def crissCrossStrat():

    angles = [[46], [50], [55], [60], [65], [70], [75], [80], [85], [89]]

    for i in range(len(angles)):

        startTime = time.time()

        boatPath = createBoatPath(coord0=[0.0, 0.0], angle0=(
            angles[i][0]*m.pi/180), stnDev=0)
        results = fishingOptimised(boatPath)

        Mean = round(sum(results)/len(results), 1)
        angles[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean, 's to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" %
              round(time.time() - startTime, 3))

    return angles


# run for different values of bS (random motion)
def randomStrat():

    stnDevs = [[0.5], [0.2], [0.1], [0.05], [0.01]]

    for i in range(len(stnDevs)):

        startTime = time.time()

        boatPath = createBoatPath(
            coord0="random", angle0="random", stnDev=stnDevs[i][0])
        results = fishingOptimised(boatPath)

        Mean = round(sum(results)/len(results), 1)
        stnDevs[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean, 's to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" %
              round(time.time() - startTime, 3))

    return stnDevs


def diamondStrat():

    tans = [[10]]

    for i in range(len(tans)):

        startTime = time.time()

        boatPath = createBoatPath(
            coord0=[0.0, L/2], angle0=m.atan(tans[i][0]), stnDev=0)
        results = fishingOptimised(boatPath)

        Mean = round(sum(results)/len(results), 1)
        tans[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean, 's to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" %
              round(time.time() - startTime, 3))

    return tans


def triangleStrat():

    tans = [[50]]

    for i in range(len(tans)):

        startTime = time.time()

        boatPath = createBoatPath(
            coord0=[0.0, 0.0], angle0=m.atan(tans[i][0]), stnDev=0)
        results = fishingOptimised(boatPath)

        Mean = round(sum(results)/len(results), 1)
        tans[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean, 's to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" %
              round(time.time() - startTime, 3))

    return tans


def spiralStrat(periods=[450], ratios=[0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]):

    # initialise
    periods = [None]+periods
    ratios = [periods]+[[Q] for Q in ratios]

    for I in range(1, len(ratios[0])):

        for i in range(1, len(ratios)):
            startTime = time.time()

            boatPath = circularMotion(constRad=False, radTimePeriod=ratios[0][I],
                                      RATIO=ratios[i][0])
            results = fishingOptimised(boatPath)

            Mean = round(sum(results)/len(results), 1)
            ratios[i].append(Mean)
            print(f"\nFrom {N} iterations,\nMean =",
                  Mean, 's to catch the fish\n')

            print("--- Computation Time: %s seconds ---\n" %
                  round(time.time() - startTime, 3))

    return ratios

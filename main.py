#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 13:46:40 2021

@author: felixdubicki-piper
"""
#MAIN#

import time
import csv
from simulation import *
from strategies import *
from extrafun import *
from config import *

startTime = time.time()

if manyBoatPaths:
    # for randomly chaging boat paths
    results = fishing()
else:
    # for same boat repeated
    boatPath = createBoatPath()
    results = fishingOptimised(boatPath)

Mean = round(sum(results)/len(results), 1)

print(f"\nFrom {N} iterations,\nMean =", Mean, 's to catch the fish\n')

print("--- Computation Time: %s seconds ---\n" %
      round(time.time() - startTime, 3))

if saveResults:
    fileName = input("Create a filename for results: ").strip()
    with open(f"{fileName}.csv", 'w') as f:
        write = csv.writer(f)

        write.writerows([i] for i in results)

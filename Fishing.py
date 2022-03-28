#!/usr/bin/env python3
# fish model
"""
Created on Sun Jan 31 13:46:40 2021

@author: felixdubicki-piper
"""

import numpy as np
import numpy.random as r
import math as m
import matplotlib.pyplot as plt
import time
import random
import csv

# Results Options ---------------------------------------------------------------

N = 1 # No. of iterations (to produce N results, higher -> more accurate)
saveResults = False
manyBoatPaths = False # New boat path every iteration (False -> faster results)



# Variables (stuff you can mess with) -------------------------------------------


A = 1000**2 # lake area, m^2 (~1000^2 m^2 average lake size, but takes long to run, min 300**2)
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



# Functions (can ignore all this) -----------------------------------------------

def startCoord(coord0):
    
    if coord0 == "random":
        coord = np.array([L*r.random(),L*r.random()])
    else:
        coord = np.array(coord0)
    
    return coord



def startAngle(angle0):
    
    if angle0 == "random":
        angle = 2*m.pi*r.random()
    else:
        angle = angle0
    
    return angle



# random & linear movement algorithm
def fish(col, fCoord, angle, t, speed, stnDev=fS, l=L, plot=plot):
    
    # generate random angle (flat distribution)
    if nrmDst == False and t>0:
        if r.random() < p: # probability of change in direction
            angle = 2*m.pi*r.random()
    
    # generate random angle (normal distribution)
    if nrmDst == True and t>0:
        angle += m.pi*r.normal(0,stnDev)
    
    if t > 0:
        # new direction
        Dir = np.array( [m.cos(angle), m.sin(angle)] ) # unit vector
        
        # new coordinate
        fCoord += speed*Dir*T
    
        # conditions at lake border (ie prevent crossing)
        if fCoord[0] <= 0 + (L-l)/2:
            fCoord[0] = 2*(L-l)/2 - fCoord[0]
            angle = m.pi-angle
        
        if fCoord[1] <= 0 + (L-l)/2:
            fCoord[1] = 2*(L-l)/2 - fCoord[1]
            angle = -angle
    
        if fCoord[0] >= (L+l)/2:
            fCoord[0] = (2*(l+L)/2 - fCoord[0])
            angle = m.pi-angle
        
        if fCoord[1] >= (L+l)/2:
            fCoord[1] = (2*(l+L)/2 - fCoord[1])
            angle = -angle
    
    # scatter plot
    '''
    if plot:
        
        plt.xlim(-1,L)
        plt.ylim(-1,L)
        
        plt.xlabel('x, $metres$')
        plt.ylabel('y, $metres$')
        
        plt.axes().set_aspect('equal')
        
        if nrmDst:
            plt.title(f'Norm Dist,  $\sigma_f$ = {fS},  $\sigma_b$ = {bS}')
        else:
            plt.title(f'Flat Dist, P = {p}')
        plt.scatter(fCoord[0],fCoord[1], s=0.1, color=col)
    '''
    
    return fCoord, angle
    

# Fishing Strategy Model 1 ------------------------------------------------------
"""
Here, the boat uses the same algorithm as the fish, so changing bS will achieve
different levels of randomness, the start angle and coords can also be changed.

For linear motion set bS = 0.

This function generates a new boat boat path each time, so be aware of how long
it takes to run.
"""

def fishing():

    results = []
    
    for n in range(N):
        
        # start time
        t = 0 # 
        
        # compute starting values
        fState = [ 'blue', startCoord(fCoord0), startAngle(fAngle0) ]
        bState = [ 'red', startCoord(bCoord0), startAngle(bAngle0) ]
        
        displacement = fState[1] - bState[1] # boat to fish vector
        distance = m.sqrt( displacement[0]**2 + displacement[1]**2 ) # boat to fish distance
    
        
        # Run simulation
        while distance > (fRad + bRad):
            
            displacement = fState[1] - bState[1] # boat to fish vector
            distance = m.sqrt( displacement[0]**2 + displacement[1]**2 ) # boat to fish distance
            
            # fish movement
            fState[1], fState[2] = fish(fState[0], fState[1], fState[2], t, f)
            
            # boat movement
            bState[1], bState[2] = fish(bState[0], bState[1], bState[2], t, b, bS, bL)
            
            '''
            # animated plot
            if animated == True and plot == True:
                plt.show() # show each step (ie animated)
                time.sleep( T*ps**-1 ) # real time delay
            
            t += T # increase time
            t = round(t,1)
            '''
            
        results.append(t)
        
        '''
        plt.legend()
        plt.show()
        '''
        
        if n%100 == 0:
            print(int(100*n/N),"%", end=", ")
    
    return results




# Fishing Strategy 2 (Linear) ---------------------------------------------------
"""
This function is more efficient as it creates a predefined boat trajectory,
rather than recalculating it for every iteration N. 
Good for linear motion (bS=0), possibly less accurate for random (bS>0).
"""

def createBoatPath(coord0=bCoord0, angle0=bAngle0, stnDev=bS,speed=b):
    
    # initialise
    t = 0
    bState = [ 'red', startCoord(coord0), startAngle(angle0) ]
    boatPath = []
    
    # set t upper bound at least 10x expected fish catching time
    while t < 500000:
        
        # new boat coord
        bState[1], bState[2] = fish(bState[0], bState[1], bState[2], 
                                    t, speed, stnDev, l=bL, plot=False)
        boatPath.append(np.array(list(bState[1])))
        
        t += T # increase time
        t = round(t,1)
        
    return boatPath



def fishingOptimised(boatPath, label1=flabel,label2=blabel,title=title):

    results = []
    
    for n in range(N):
        
        # start time
        t = 0
        
        # initialise record of fish and boat paths (for plotting)
        fishAndBoatPath = [[],[]]
        
        # compute starting values
        fState = [ 'blue', startCoord(fCoord0), startAngle(fAngle0) ]
        
        # inital boat to fish vector
        displacement = fState[1] - boatPath[0] 
        # initial boat to fish distance
        distance = m.sqrt( displacement[0]**2 + displacement[1]**2 ) 
    
        
        # Run simulation
        for i in range(len(boatPath)):
            
            # new boat to fish vector
            displacement = fState[1] - boatPath[i]
            # new boat to fish distance
            distance = m.sqrt( displacement[0]**2 + displacement[1]**2 )
            
            # fish movement
            fState[1], fState[2] = fish(fState[0], fState[1], fState[2], t, f)
            
            if plot:
                # save current coord. of fish and boat
                fishAndBoatPath[0].append(np.array( list(fState[1]) ))
                fishAndBoatPath[1].append(boatPath[i])
            
            # animated plot
            if animated:
                plt.xlim(-1,L)
                plt.ylim(-1,L)
                
                plt.xlabel('x, $metres$')
                plt.ylabel('y, $metres$')
                
                plt.axes().set_aspect('equal')
                
                plt.title(title)
                
                if scatter:
                    plt.scatter(fState[1][0],fState[1][1], s=0.1, color='b')
                    plt.scatter(boatPath[i][0], boatPath[i][1], s=0.1, color='r')
                else:
                    fishX = [ coord[0] for coord in fishAndBoatPath[0] ]
                    fishY = [ coord[1] for coord in fishAndBoatPath[0] ]
                    
                    boatX = [ coord[0] for coord in fishAndBoatPath[1] ]
                    boatY = [ coord[1] for coord in fishAndBoatPath[1] ]
                    
                    plt.plot(fishX, fishY, c='b', label='$\sigma_{fish}$ = '+str(fS))
                    plt.plot(boatX, boatY, c='r', label='$\sigma_{boat}$ = '+str(bS))
                    # plt.plot(fishAndBoatPath[0][i][0],fishAndBoatPath[0][i][1], color='b')
                    # plt.plot(fishAndBoatPath[1][i][0],fishAndBoatPath[1][i][1], color='r')
                
                plt.show() # show each step (ie animated)
                # time.sleep( T*ps**-1 ) # real time delay
            
            
            if distance <= (fRad + bRad):
                break
            
            t += T # increase time
            t = round(t,1)
        
        results.append(t)
        
        
        # plot fish and boat paths
        if plot:
            fishX = [ coord[0] for coord in fishAndBoatPath[0] ]
            fishY = [ coord[1] for coord in fishAndBoatPath[0] ]
            
            boatX = [ coord[0] for coord in fishAndBoatPath[1] ]
            boatY = [ coord[1] for coord in fishAndBoatPath[1] ]
            
            fig = plt.figure(dpi=300, figsize=(2.5, 2.5))
            ax = fig.add_subplot(111)
            
            if scatter:
                ax.scatter(fishX, fishY, s=0.1, c='b', label='$\sigma_{fish}$ = '+str(fS))
                ax.scatter(boatX, boatY, s=0.1, c='r', label='$\sigma_{boat}$ = '+str(bS))
            
            else:
                ax.plot(fishX, fishY, c='b', label='$\sigma_{fish}$ = '+str(fS))
                ax.plot(boatX, boatY, c='r', label='$\sigma_{boat}$ = '+str(bS))
            
            ax.set_aspect('equal','box')
            ax.set(xlim = (-1,L), ylim = (-1,L) )
            
            # plt.legend(bbox_to_anchor=(0.01, 0.905, 0.98, 0.1), loc='lower left',
            #            ncol=2, mode="expand", borderaxespad=0., markerscale=10, 
            #            framealpha=1)
            
            plt.legend(loc='upper right', markerscale=10, fontsize=8, framealpha=0.95,
                       handletextpad=0.2, handlelength=1.)
            
            ax.set_xlabel('X, $metres$', fontsize=8)
            ax.set_ylabel('Y, $metres$', fontsize=8)
            
            ax.tick_params(axis='x', which='major', labelsize=8)
            ax.tick_params(axis='y', which='major', labelsize=8, labelrotation=90)
            
            ax.set_title(title, y=1.0)
        
            plt.show()
        
        if n>0 and n%100 == 0:
            print(int(100*n/N),"%", end=", ")
        
    return results



# plot Precreated boat path
def plotBoat(boatPath, tLim=1800,col='r',title=title,label='Boat'):
    
    indexLim = round(tLim/T)
    
    boatX = [coord[0] for coord in boatPath]
    boatY = [coord[1] for coord in boatPath]
    
    fig = plt.figure(dpi=300, figsize=(2.5, 2.5))
    ax = fig.add_subplot(111)
    
    if scatter:
        ax.scatter(boatX[:indexLim], boatY[:indexLim], s=0.1, c=col, label=label)
    else:
        ax.plot(boatX[:indexLim], boatY[:indexLim], c=col, label=label)
    
    
    ax.set_aspect('equal','box')
    ax.set(xlim = (0,L), ylim = (0,L) )
    
    plt.legend(loc='upper right', markerscale=10, fontsize=8, framealpha=0.95,
               handletextpad=0.2, handlelength=1.)
    
    ax.set_xlabel('X, $metres$', fontsize=8)
    ax.set_ylabel('Y, $metres$', fontsize=8)
    
    ax.tick_params(axis='x', which='major', labelsize=8)
    ax.tick_params(axis='y', which='major', labelsize=8, labelrotation=90)
    
    ax.set_title(title, y=1.0)
        
    plt.show()



# run for different starting angles (betwen 45º and 90º exclusive)
def crissCrossStrat():
    
    angles = [ [46], [50], [55], [60], [65], [70], [75], [80], [85], [89] ]
    
    for i in range(len(angles)):
        
        startTime = time.time()
        
        boatPath = createBoatPath( coord0=[0.0,0.0], angle0=(angles[i][0]*m.pi/180), stnDev=0 )
        results = fishingOptimised(boatPath)
        
        Mean = round(sum(results)/len(results),1)
        angles[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean,'s to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" % round(time.time() - startTime,3) )
        
    return angles



# run for different values of bS (random motion)
def randomStrat():
    
    stnDevs = [ [0.5], [0.2], [0.1], [0.05], [0.01] ]
    
    for i in range(len(stnDevs)):
        
        startTime = time.time()
        
        boatPath = createBoatPath( coord0="random", angle0="random", stnDev=stnDevs[i][0] )
        results = fishingOptimised(boatPath)
        
        Mean = round(sum(results)/len(results),1)
        stnDevs[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean,'s to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" % round(time.time() - startTime,3) )
        
    return stnDevs


def diamondStrat():
    
    tans = [ [10] ]
    
    for i in range(len(tans)):
        
        startTime = time.time()
        
        boatPath = createBoatPath( coord0=[0.0,L/2], angle0 = m.atan( tans[i][0] ), stnDev=0 )
        results = fishingOptimised(boatPath)
        
        Mean = round(sum(results)/len(results),1)
        tans[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean,'s to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" % round(time.time() - startTime,3) )
        
    return tans

def triangleStrat():
    
    tans = [ [50] ]
    
    for i in range(len(tans)):
        
        startTime = time.time()
        
        boatPath = createBoatPath( coord0=[0.0,0.0], angle0 = m.atan( tans[i][0] ), stnDev=0 )
        results = fishingOptimised(boatPath)
        
        Mean = round(sum(results)/len(results),1)
        tans[i].append(Mean)
        print(f"\nFrom {N} iterations,\nMean =", Mean,'s to catch the fish\n')

        print("--- Computation Time: %s seconds ---\n" % round(time.time() - startTime,3) )
        
    return tans



def circularMotion(constRad=True, RATIO=0.95, radTimePeriod=450, speed=b,animated=False):
    
    # initialise
    t = 0
    bCoord = startCoord(bCoord0)
    boatPath = []
    radRateChange = 2*(m.pi)/(2*radTimePeriod)
    invRATIO = 1-RATIO
    jumpStep = m.pi/(2*radRateChange)
    
    # define start/max radius
    RAD = m.sqrt( (bCoord[0]-L/2)**2+(bCoord[1]-L/2)**2 )
    
    # set t upper bound at least 10x expected fish catching time
    while t < 500000:
        
        if constRad:
            bCoord = [ RAD*m.cos(t*speed/RAD) + L/2, RAD*m.sin(t*speed/RAD) + L/2]
        else:
            # determine current radius
            rad = RAD*( RATIO*m.cos(radRateChange*t)**2 + 1-RATIO)
            # determine number of angle steps
            nStep = (1+t//jumpStep)//2
            
            # determine new angle
            angleFunc = (speed/RAD)*1/(radRateChange*m.sqrt(invRATIO))*m.atan( m.sqrt(invRATIO)*m.tan(radRateChange*t) )
            # angle step, as function 'resets' periodcally, causing jumps in the trajectory
            angleFunc += nStep*(speed/RAD)*1/(radRateChange*m.sqrt(invRATIO))*m.pi
            
            # new boat coord
            bCoord = [ rad*m.cos(angleFunc) + L/2, rad*m.sin(angleFunc) + L/2]
        
        boatPath.append(bCoord)
        
        if animated:
        
            plt.xlim(-1,L)
            plt.ylim(-1,L)
            
            plt.xlabel('x, $metres$')
            plt.ylabel('y, $metres$')
            
            plt.axes().set_aspect('equal')
            
            plt.title('Spiral')
            
            plt.scatter(bCoord[0],bCoord[1], s=0.1, color='red')
            # time.sleep( T*ps**-1 ) # real time delay
            
            plt.show()
        
        t += T # increase time
        t = round(t,1)
        
    return boatPath
    


def spiralStrat(periods=[450], ratios=[0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4] ): 

    # initialise
    periods = [None]+periods
    ratios = [periods]+[[Q] for Q in ratios]

    for I in range( 1, len(ratios[0]) ):
        
        for i in range( 1,len(ratios) ):
            startTime = time.time()
                
            boatPath = circularMotion(constRad=False, radTimePeriod=ratios[0][I], 
                                      RATIO=ratios[i][0])
            results = fishingOptimised(boatPath)
            
            Mean = round(sum(results)/len(results),1)
            ratios[i].append(Mean)
            print(f"\nFrom {N} iterations,\nMean =", Mean,'s to catch the fish\n')
        
            print("--- Computation Time: %s seconds ---\n" % round(time.time() - startTime,3) )
    
    return ratios



# %% MAIN -----------------------------------------------------------------------

startTime = time.time()

if manyBoatPaths:
    results = fishing()
else:
    boatPath = createBoatPath()
    results = fishingOptimised(boatPath)
    
Mean = round(sum(results)/len(results),1)

print(f"\nFrom {N} iterations,\nMean =", Mean,'s to catch the fish\n')

print("--- Computation Time: %s seconds ---\n" % round(time.time() - startTime,3) )

if saveResults:
    fileName = input("Create a filename for results: ").strip()
    with open(f"{fileName}.csv", 'w') as f:
        write = csv.writer(f)
        
        write.writerows([i] for i in results)
        
        
# %% Additional Mini Functions --------------------------------------------------

def plotFishTrajectory(tLim=18000, plotNum=1,σ=fS,speed=f):
    for i in range(plotNum):
        fishPath=createBoatPath(coord0=fCoord0,angle0=fAngle0,stnDev=σ, speed=speed)
        plotBoat(fishPath, tLim=tLim, col='b',label=f"Fish\nt={tLim}s", title="Fish Trajectory, "+f"$\sigma={σ},$ $f={speed}$ $m/s$")
        
def plotBoatSpiral(tLim=1800,Q=0.6,period=450):
    boatPath = circularMotion(constRad=False, RATIO=(1-Q), radTimePeriod=period, speed=b)
    plotBoat(boatPath, tLim=tLim, col='r',label=f"Boat\nt={tLim}s", title="Spiral, $T_{period}=$"+str(period)+"s, $Ratio=$"+str(Q))
    
def plotBoatTrajectory(tLim=1800, xy0=[0.0,0.0], α0=m.pi/3, σ=0,title=title, label=blabel ):
    boatPath=createBoatPath(coord0=xy0,angle0=α0,stnDev=σ, speed=b)
    plotBoat(boatPath, tLim=tLim, col='r',label=label, title=title)


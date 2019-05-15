#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
from subprocess import Popen, list2cmdline, PIPE
import numpy as np
import os.path
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from matplotlib import cm

# this is a helper script to plot results (if needed)
# just provided but some parameters need to be changed

# nx should be equal to Nx-1 where Nx is the value set in main.cc
# this is becuase of the convention we use to discretize
nx = 64
ny = 31;

startFrom = int(sys.argv[1])
endAt = int(sys.argv[2])
freq = int(sys.argv[3])

def getXY():
    dd = np.loadtxt("xy.txt")
    x,y = dd[:,0], dd[:,1]
    x = x.reshape(ny,nx)
    y = y.reshape(ny,nx)
    return x,y

x,y = getXY()
for i in range(startFrom, endAt, freq):
    d0 = np.loadtxt("sol_"+str(i)+".txt")
    c0 = d0[0:-1:4]  #T
    c1 = d0[1:-1:4]  #H2
    c2 = d0[2:-1:4]  #O2
    #c3 = d0[3:len(d0)+1:4]

    c0 = c0.reshape(ny,nx)
    c1 = c1.reshape(ny,nx)
    c2 = c2.reshape(ny,nx)

    nLevs = 25

    fig = plt.figure(1)
    #ax1 = fig.add_subplot(141)
    #ax2 = fig.add_subplot(142)
    #ax3 = fig.add_subplot(143)

    #lev1 = np.linspace(np.min(c0), np.max(c0), nLevs)
    plt.contourf(x,y,c0, cmap=cm.jet)
    #plt.set_aspect(aspect=1)

    # #lev2 = np.linspace(np.min(c1), np.max(c1), nLevs)
    #ax2.contourf(x,y,c1, cmap=cm.jet)
    #ax2.set_aspect(aspect=1)

    # #lev3 = np.linspace(np.min(c2), np.max(c2), nLevs)
    # ax3.contourf(x,y,c2, cmap=cm.jet)
    # ax3.set_aspect(aspect=1)

    plt.clim(0,1600)
    plt.colorbar()
    plt.pause(0.01)
    plt.show()

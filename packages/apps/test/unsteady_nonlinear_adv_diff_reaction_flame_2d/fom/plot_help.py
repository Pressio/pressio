#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
from subprocess import Popen, list2cmdline, PIPE
import numpy as np
import os.path
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from matplotlib import cm

# nx should be equal to Nx-1 where Nx is the value set in main.cc
# this is becuase of the convention we use to discretize
nx = 72
ny = 36

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
    c3 = d0[3:len(d0)+1:4] #H2O

    c0 = c0.reshape(ny,nx)
    c1 = c1.reshape(ny,nx)
    c2 = c2.reshape(ny,nx)
    c3 = c3.reshape(ny,nx)

    nLevs = 35
    fig = plt.figure(1)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    lev1 = np.linspace(np.min(c0), np.max(c0), nLevs)
    ax1.imshow(c0, cmap=cm.jet, origin='lower')#,interpolation='bicubic')
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    #plt.contourf(x,y,c0, lev1,cmap=cm.jet)
    #plt.plot(x,y,'ok', markersize=5)
    #ax = plt.gca()
    ax1.set_aspect(aspect=1)

    ax2.imshow(c1, cmap=cm.brg, origin='lower')#,interpolation='bicubic')
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)

    ax3.imshow(c2, cmap=cm.brg, origin='lower')#,interpolation='bicubic')
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)

    ax4.imshow(c3, cmap=cm.brg, origin='lower')#,interpolation='bicubic')
    ax4.get_xaxis().set_visible(False)
    ax4.get_yaxis().set_visible(False)

    #plt.clim(300,1800)
    # plt.xlim(0,1.8)
    # plt.ylim(0,0.9)
    #plt.colorbar()
    plt.pause(0.01)
    plt.show()

#! /bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os  
import sys

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


searchdir = sys.argv[1]#'8x8/Ref'
#searchdir = '80x80-03-09-15--10-03'
#searchdir = '80x80-03-05-15--20-28'
#searchdir = 'comm_stats_6400'
nrows = int(sys.argv[2])
ncols = int(sys.argv[3])
cmapname = 'nipy_spectral'

x = np.array(range(1,nrows,1))
y = np.array(range(1,ncols,1))

grid = np.zeros((nrows,ncols))


fig = plt.figure(1)
grid = np.zeros((nrows,ncols))
plt.title(searchdir)

first = 0
with open(searchdir) as f:
    for line in f:
        if first == 0:
            receiver,sizei = line.split()
            xi = int(receiver)%ncols
            yi = int(receiver)/ncols
            sizei = float(sizei)
            grid[xi,yi] += sizei
        first = 0

plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
           interpolation='nearest', cmap=plt.get_cmap(cmapname))
plt.colorbar()
 

fig.savefig(searchdir + ".pdf", format='pdf')

#### plt.title('Total comm. volume')
####
####first = 0
####with open("%s/total.dat" % searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
####plt.colorbar()
#### 
####
####fig.savefig("total" + ".pdf", format='pdf')
####
####
####fig = plt.figure(2)
####plt.title('Avg. comm. volume')
####
####first = 0
####grid = np.zeros((nrows,ncols))
####with open("%s/avg.dat" % searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
####plt.colorbar()
####
####
####fig.savefig("avg" + ".pdf", format='pdf')
####
####
####
####fig = plt.figure(3)
####grid = np.zeros((nrows,ncols))
####plt.title('Total comm. volume Reduce L')
####
####first = 0
####with open("%s/reduceL_total.dat"%searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
####plt.colorbar()
#### 
####
####fig.savefig("reduceL_total" + ".pdf", format='pdf')
####
####
#####fig = plt.figure(4)
#####grid = np.zeros((nrows,ncols))
#####plt.title('Total comm. volume send L CD')
#####
#####first = 0
#####with open("%s/send_L_CD_total.dat"%searchdir) as f:
#####    for line in f:
#####        if first == 0:
#####            receiver,sizei = line.split()
#####            xi = int(receiver)%ncols
#####            yi = int(receiver)/ncols
#####            sizei = float(sizei)
#####            grid[xi,yi] += sizei
#####        first = 0
#####
#####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
#####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
#####plt.colorbar()
##### 
#####
#####fig.savefig("send_L_CD_total" + ".pdf", format='pdf')
#####
#####fig = plt.figure(5)
#####grid = np.zeros((nrows,ncols))
#####plt.title('Total comm. volume bcast L')
#####
#####first = 0
#####with open('comm_stats_6400/bcastL_total.dat') as f:
#####    for line in f:
#####        if first == 0:
#####            receiver,sizei = line.split()
#####            xi = int(receiver)%ncols
#####            yi = int(receiver)/ncols
#####            sizei = float(sizei)
#####            grid[xi,yi] += sizei
#####        first = 0
#####
#####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
#####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
#####plt.colorbar()
##### 
#####
#####fig.savefig("bcastL_total" + ".pdf", format='pdf')
####
####
####fig = plt.figure(6)
####grid = np.zeros((nrows,ncols))
####plt.title('Total comm. volume bcast U')
####
####first = 0
####with open("%s/bcastU_total.dat"%searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
####plt.colorbar()
#### 
####
####fig.savefig("bcastU_total" + ".pdf", format='pdf')
####
####fig = plt.figure(7)
####grid = np.zeros((nrows,ncols))
####plt.title('Total comm. volume reduce L (sender)')
####
####first = 0
####with open("%s/sender_reduceL_total.dat"%searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
####plt.colorbar()
#### 
####
####fig.savefig("sender_reduceL_total" + ".pdf", format='pdf')
####
####
####
####fig = plt.figure(8)
####grid = np.zeros((nrows,ncols))
####plt.title('Total comm. volume bcast U (sender)')
####
####first = 0
####with open("%s/sender_bcastU_total.dat"%searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####print np.amin(grid)
####print np.amax(grid)
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname))
####plt.colorbar()
#### 
####
####fig.savefig("sender_bcastU_total" + ".pdf", format='pdf')
####
####fig = plt.figure(9)
####grid = np.zeros((nrows,ncols))
####plt.title('Total comm. volume (sender)')
####
####first = 0
####with open("%s/sender_total.dat"%searchdir) as f:
####    for line in f:
####        if first == 0:
####            receiver,sizei = line.split()
####            xi = int(receiver)%ncols
####            yi = int(receiver)/ncols
####            sizei = float(sizei)
####            grid[xi,yi] += sizei
####        first = 0
####
####
####print np.amin(grid)
####print np.amax(grid)
####plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
####           interpolation='nearest', cmap=plt.get_cmap(cmapname),vmin=np.amin(grid), vmax=np.amax(grid))
####plt.colorbar()
#### 
####
####fig.savefig("sender_total" + ".pdf", format='pdf')
####
####
#####with open("%s/sender_bcastU_msg.dat"%searchdir) as f:
#####    fig = plt.figure(10)
#####    grid = np.zeros((nrows,ncols))
#####    plt.title('Total number of messages (sender)')
#####    first = 0
#####    for line in f:
#####        if first == 0:
#####            receiver,sizei = line.split()
#####            xi = int(receiver)%ncols
#####            yi = int(receiver)/ncols
#####            sizei = float(sizei)
#####            grid[xi,yi] += sizei
#####        first = 0
#####
#####
#####    print np.amin(grid)
#####    print np.amax(grid)
#####    plt.imshow(grid, extent=(x.min()-0.5, x.max()+0.5, y.max()+0.5, y.min()-0.5),
#####               interpolation='nearest', cmap=plt.get_cmap(cmapname),vmin=np.amin(grid), vmax=np.amax(grid))
#####    plt.colorbar()
#####     
#####    
#####    fig.savefig("sender_bcastU_msg" + ".pdf", format='pdf')
####
####
####

plt.show()

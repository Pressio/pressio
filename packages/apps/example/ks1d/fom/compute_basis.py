import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl

dt = 0.1
t0 = 0.0
L = 128.0

ss_rate = 1 

snapshots = np.loadtxt("primal.dat")

nnode = snapshots.shape[0]
Nsamp = snapshots.shape[1]


# sample as needed:
snapshots = snapshots[:,0::ss_rate]

# subtract IC
ic= snapshots[:,0]
snapshots = snapshots[:,1:] - ic[:,np.newaxis] 


# compute SVD
U,sigma,_ = spl.svd(snapshots,full_matrices=False)

energy = np.cumsum(sigma**2)/np.sum(sigma**2)


# save to file
np.savetxt("restart.inp",ic)
np.savetxt("basis.txt",U)

plt.figure()

plt.plot(energy)

plt.show()


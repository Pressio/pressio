
import numpy as np
import matplotlib.pyplot as plt

numPartitions = 8
numRanks = 5

fomTotalDofsCount = 0
rows_indices = {}
bases = {}
for i in range(numPartitions):
  bases[i] = np.loadtxt("bases_p_" + str(i) + ".txt", dtype=float)
  rows = np.loadtxt("state_vec_rows_wrt_full_mesh_p_" + str(i) + ".txt", dtype=int)
  fomTotalDofsCount += len(rows)
  rows_indices[i] = rows.copy()
print(fomTotalDofsCount)

romY = np.array([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8])

phi = np.zeros((fomTotalDofsCount, len(romY)))
for p in range(numPartitions):
  A = bases[p]
  rows = rows_indices[p]
  for k in range(len(rows)):
    phi[rows[k], p*3:p*3+3] = A[k, :].copy()

fomY = np.dot(phi, romY)
np.savetxt("py_fom_state.txt", fomY, fmt='%14.12f')

fig = plt.figure()
plt.plot(np.arange(fomTotalDofsCount), fomY, '*k', markersize=5)

for i in range(numRanks):
  y = np.loadtxt("fomTmpState_"+str(i)+".txt", dtype=float)
  gids = np.loadtxt("tpetra_map_gids_"+str(i)+".txt", dtype=int)
  plt.plot(gids, y, 'ob', markerfacecolor='None', markersize=6)
plt.show()

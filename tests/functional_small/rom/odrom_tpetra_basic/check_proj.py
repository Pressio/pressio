
import numpy as np

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

# assemble global basis
phi = np.zeros((fomTotalDofsCount, numPartitions*3))
for p in range(numPartitions):
  A = bases[p]
  rows = rows_indices[p]
  for k in range(len(rows)):
    phi[rows[k], p*3:p*3+3] = A[k, :].copy()


velo_per_rank = {}
tp_gids_per_rank = {}
for r in range(5):
  velo_per_rank[r] = np.loadtxt("velo_"+str(r)+".txt", dtype=float)
  tp_gids_per_rank[r] = np.loadtxt("tpetra_map_gids_"+str(r)+".txt", dtype=int)

fomVelo = np.zeros(fomTotalDofsCount)
for r in range(5):
  rankVelo = velo_per_rank[r]
  for k in range(len(rankVelo)):
    gid = tp_gids_per_rank[r][k]
    fomVelo[gid] = rankVelo[k]

y = np.dot(phi.transpose(), fomVelo)
print(y)

for i in [0,1,2,3,4]:
  computedRomState = np.loadtxt("rom_state_projected_"+str(i)+".txt")
  assert(np.allclose(computedRomState, y,rtol=1e-8, atol=1e-10))


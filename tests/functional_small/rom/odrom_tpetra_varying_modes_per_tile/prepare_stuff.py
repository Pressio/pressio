
import numpy as np

# ============================================
# for each rank, need to fill the cell gids
# ============================================
def _mapCellGidsToStateDofsGids(pCellGids):
  myl = [None]*len(pCellGids)*4
  myl[0::4] = [int(i*4)   for i in pCellGids]
  myl[1::4] = [int(i*4+1) for i in pCellGids]
  myl[2::4] = [int(i*4+2) for i in pCellGids]
  myl[3::4] = [int(i*4+3) for i in pCellGids]
  return np.array(myl)

cellGids = {}

cellGids[0] = [*range(0, 12, 1)]
for i in range(1,23):
  cellGids[0] += [*range(i*46, i*46+12, 1)]
print(cellGids[0])
np.savetxt("cell_gids_on_rank_0.txt", cellGids[0], fmt='%d')
np.savetxt("tpetra_map_gids_0.txt", _mapCellGidsToStateDofsGids(cellGids[0]), fmt='%d')
print("")

cellGids[1] = [*range(380, 395, 1)]
for i in range(1,15):
  startat = i*46+380
  cellGids[1] += [*range(startat, startat+15, 1)]
print(cellGids[1])
np.savetxt("cell_gids_on_rank_1.txt", cellGids[1], fmt='%d')
np.savetxt("tpetra_map_gids_1.txt", _mapCellGidsToStateDofsGids(cellGids[1]), fmt='%d')
print("")

cellGids[2] = [*range(12, 35, 1)]
for i in range(1,8):
  startat = i*46+12
  cellGids[2] += [*range(startat, startat+23, 1)]
print(cellGids[2])
np.savetxt("cell_gids_on_rank_2.txt", cellGids[2], fmt='%d')
np.savetxt("tpetra_map_gids_2.txt", _mapCellGidsToStateDofsGids(cellGids[2]), fmt='%d')
print("")

cellGids[3] = [*range(395, 414, 1)]
for i in range(1,15):
  startat = i*46+395
  cellGids[3] += [*range(startat, startat+19, 1)]
print(cellGids[3])
np.savetxt("cell_gids_on_rank_3.txt", cellGids[3], fmt='%d')
np.savetxt("tpetra_map_gids_3.txt", _mapCellGidsToStateDofsGids(cellGids[3]), fmt='%d')
print("")

cellGids[4] = [*range(35, 46, 1)]
for i in range(1,8):
  startat = i*46+35
  cellGids[4] += [*range(startat, startat+11, 1)]
print(cellGids[4])
np.savetxt("cell_gids_on_rank_4.txt", cellGids[4], fmt='%d')
np.savetxt("tpetra_map_gids_4.txt", _mapCellGidsToStateDofsGids(cellGids[4]), fmt='%d')
print("")

# ============================================
# fake bases
# ============================================
numTiles = 8
basesPerTile = [3,2,4,5,2,3,4,3]
np.savetxt("num_modes_per_tile.txt", basesPerTile, fmt='%3d')

for i in range(numTiles):
  filen = "state_vec_rows_wrt_full_mesh_p_" + str(i) + ".txt"
  numStencilDofs = len(np.loadtxt(filen))
  print(numStencilDofs)
  M = np.random.random((numStencilDofs, basesPerTile[i]))
  np.savetxt("bases_p_"+str(i)+".txt", M)

totModesCount = np.sum(basesPerTile)
np.savetxt("rom_state_rand.txt", np.random.random(totModesCount))

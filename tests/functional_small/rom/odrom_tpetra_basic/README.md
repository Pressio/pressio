
## "rank_partition_mapping.txt"

  - info about which tiles are owned by which rank, example
  ```
  5
  0 0  4 -1 -1
  1 1  2  5  6
  2 1  2 -1 -1
  3 2  3  6  7
  4 3 -1 -1 -1
  ```
  where
   - 5 on first line specifies how many ranks exist
   - each row contains the rank (on first column) and list of partitions ids belonging to that rank
   so
   ```
     0 0  4 -1 -1
   ```
   is for rank==0, which "owns" partitions 0 and 4
   and
   ```
   1  1  2  5  6
   ```
   this is the line for rank=1 and we say that tiles 1,2,5,6 are partially or fully owned by rank=1

   - note that for some rows you have to pad with -1 to ensure all rows have same # of values


## bases_p_{0,1,...,n}.txt

  - contains bases in ascii format for partition `0` to `n`:

  - if i-th partitions has `N_i` state dofs, then the file should have a matrix of `N_i \times K` where `K` is num of modes

  - example for `N_i = 12` (for example a partition with 3 cells and 4 dofs/cell) and `K=3`:
    ```
    7.385650552745359754e-01 9.558753197446867578e-02 1.287464274747743831e-01
    1.364995925095853213e-01 5.684999169722433354e-01 9.270582132695502908e-01
    1.149458688516840077e-01 8.244847250788629456e-01 6.702373153885047286e-01
    3.420770700693368527e-01 7.156607043284927139e-01 4.942478840988473454e-01
    3.597777058002636918e-01 4.069035497079750430e-01 5.987705696120064758e-01
    2.847022641872032356e-01 4.844329677175417403e-01 6.249961995866215592e-01
    2.140431295489493291e-01 6.308919328646183100e-01 1.392343353778845438e-01
    1.913842174742834690e-01 3.005851445767541152e-01 1.353245236166178422e-01
    3.191367339286875771e-01 3.462379506471334745e-01 2.451953765270737939e-01
    2.063414187655948639e-01 8.687565068637496113e-01 3.652977291299257523e-01
    4.703000961625675158e-01 6.570333646605375222e-01 3.436119552021041912e-01
    8.544088020812895534e-01 6.657001343987848374e-01 3.291331628879936577e-01
    ```

## state_vec_rows_wrt_full_mesh_p_{0,1,...,n}.txt

  - contains the GLOBAL gids of the FOM state vector that belong to partition `0` to `n`
  - for each partition, this allows us to know which entries in the FOM state vector belong to it

  - example:
  ```
  30
  31
  32
  33
  34
  35
  36
  37
  38
  39
  40
  41
  42
  43
  44
  45
  46
  47
 184
 185
 186
 187
 188
 189
 190
 191
 192
 ```

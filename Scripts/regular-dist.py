#!/usr/bin/python3
import time
import numpy as np

######################################################################
## Generate nearest neighbor list from coords
## Assumes only C present
######################################################################
def find_bonds(coords,cutoff,abc):
    bonds = []
    lengths = []
    for i in range(0,len(coords)):
        ix = coords[i][0]
        iy = coords[i][1]
        iz = coords[i][2]
        NeighborElements = []
        for j in range(0,len(coords)):
            jx = coords[j][0]
            jy = coords[j][1]
            jz = coords[j][2]
            dx = abs(ix - jx)
            if dx >= (abc[0]*0.5):
                dx -= abc[0]
            dy = abs(iy - jy)
            if dy >= (abc[1]*0.5):
                dy -= abc[1]
            dz = abs(iz - jz)
            if dz >= (abc[2]*0.5):
                dz -= abc[2]
            d = (dx**2 + dy**2 + dz**2 )**0.5
            if d <= cutoff and i != j:
                bonds.extend(sorted([i,j]))
                lengths.extend([d])

    bond_topo, b_ind = np.unique(np.array(bonds), axis=0, return_index=True)
    bond_lengths = np.array(lengths)[b_ind]
    return(bond_topo,bond_lengths)


if __name__=="__main__":
    # Inputs
    natoms = 5000
    np.random.seed(666)
    positions = np.random.rand(natoms, 3)*10.0
    box = np.full(3, 10)
    pbc = np.full(3, True)
    r_cut = 1.5

    # Run the function
    tic = time.perf_counter()
    bond_topo, bond_lengths = find_bonds(positions, r_cut, box)
    toc = time.perf_counter()

    # print the results
    print(f'{natoms} atoms')
    print(f'There are {bond_topo.shape[0]} bonds')
    print(f'Bond Topology: {bond_topo}')
    print(f'Bond Lengths: {bond_lengths}')
    print(f'Result calculated in {toc - tic:0.4f} seconds')

    np.savetxt('test.xyz', positions, fmt='%5f', header=f'{natoms}\n{box[0]} {box[1]} {box[2]}', comments= '')

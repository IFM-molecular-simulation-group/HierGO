#!/usr/bin/python3
from argparse import ArgumentParser
import sys, math, os, statistics, time, random
import numpy as np


######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser(description='')
    ## --infile filename
    parser.add_argument('--infile', help='<Required> Set flag; input file \
    name (pdb)', required=True, type=str)

    ## --outfile filename
    parser.add_argument('--outfile', nargs='?', default="trimmed.pdb", \
    type=str, help="output file name (pdb)")

    ## --atoms 10
    parser.add_argument('--atoms', nargs='?', type=int, help='optionally specify \
    an additional number of atoms to remove')

    ## --orthogonalize
    parser.add_argument('--orthogonalize', action='store_true', help='convert from \
    rhombohedral cell to orthogonal cell')

    args = parser.parse_args()
    return(args.infile,args.outfile,args.atoms,args.orthogonalize)


######################################################################
## Read pdb to local structure format
######################################################################
def read_pdb_to_structure(InputFile):
   coords = []
   atomtypes = []
   x_min = 1000
   y_min = 1000
   for line in open(InputFile):
       if "CRYST1" in line:
           PBC = line.split()
           abc = [ float(PBC[1]), float(PBC[2]), float(PBC[3]), float(PBC[4]), float(PBC[5]), float(PBC[6])]
       elif "ATOM" in line or "HETATM" in line:
           fields = line.split()
           coords.append( [ float(fields[5]), float(fields[6]), float(fields[7]) ] )
           atomtypes.append( fields[2] )
           if abs(float(fields[5])) < x_min and abs(float(fields[6])) < y_min:
               originatom = int(fields[1])
               x_min = abs(float(fields[5]))
               y_min = abs(float(fields[6]))
   return(coords,atomtypes,abc,originatom)

######################################################################
## There used to be a problem with overlapping atoms before I added 
## the origin movement, I will keep this function here in case it 
## becomes a problem again
######################################################################
def remove_overlapping_atoms(coords,atomtypes,abc):
    NNList = []
    newcoords = []
    newatomtypes = []
    RemoveAtoms = []
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
            if d <= 0.01 and i != j:
                RemoveAtoms.append(j)
    half = len(RemoveAtoms)//2
    #if half > 0:
        #print("Removing %s atoms"%(half))
    for i in range(0,len(coords)):
        if i not in RemoveAtoms[:half]:
            newcoords.append(coords[i])
            newatomtypes.append(atomtypes[i])
    return(newcoords,newatomtypes)

######################################################################
## nearest neighbors for origin atom, ignores pbc for the purpose of 
## always having the origin atom in the same way
######################################################################
def single_nearest_neighbors(originatom,coords):
    cutoff = 1.7
    ix = coords[originatom][0]
    iy = coords[originatom][1]
    iz = coords[originatom][2]
    NeighborElements = []
    for j in range(0,len(coords)):
        jx = coords[j][0]
        jy = coords[j][1]
        jz = coords[j][2]
        dx = abs(ix - jx)
        dy = abs(iy - jy)
        dz = abs(iz - jz)
        d = (dx**2 + dy**2 + dz**2 )**0.5
        if d <= cutoff and originatom != j:
            NeighborElements.append(j)
    if len(NeighborElements) == 1:
        originatom = NeighborElements[0]
    moveby = [ coords[originatom][0], coords[originatom][1], coords[originatom][2] ]
    return(moveby)

######################################################################
## Determine new y length
######################################################################
def orthogonalize_graphene(PBC,Coords):
    print("Old PBC: %s"%(PBC))
    if PBC[5] != 90.0:
        #PBC[0] = PBC[0] * ( np.sin(np.deg2rad(PBC[5])) / np.sin(np.deg2rad(PBC[3])))
        PBC[1] = PBC[1] * ( np.sin(np.deg2rad(PBC[5])) / np.sin(np.deg2rad(PBC[4])))
        PBC[5] = 90
    print("New PBC: %s"%(PBC))
    for i in range(0,len(Coords)):
        if Coords[i][1] >= PBC[1]:
            Coords[i][1] -= PBC[1]
    return(PBC,Coords)

######################################################################
## Wrap coordinates back to within the periodic boundary conditions
## (ignore z)
######################################################################
def wrap_orthogonal_coords(coords,abc):
    for atom in range(0,len(coords)):
        if coords[atom][0] < 0:
            coords[atom][0] += abc[0]
        elif coords[atom][0] > abc[0]:
            coords[atom][0] -= abc[0]
        if coords[atom][1] < 0:
            coords[atom][1] += abc[1]
        elif coords[atom][1] > abc[1]:
            coords[atom][1] -= abc[1]
    return(coords)

######################################################################
## 
######################################################################
def make_consistant_origin(coords,originatom):
    moveby = single_nearest_neighbors(originatom,coords)
    for i in range(0,len(coords)):
        coords[i][0] += moveby[0]
        coords[i][1] += moveby[1]
        coords[i][2] += moveby[2]
    return(coords)


######################################################################
## FV NN
######################################################################
def distance(x0, x1, box, pbc):
    # xo is a position of one atom, x1 is an array of positions
    # use the pbc bool mask to set the periodicity
    delta = np.abs(x0 - x1)
    delta[:,pbc] -= box[pbc] * np.round(delta[:,pbc]/(box[pbc]))
    return np.sqrt((delta ** 2).sum(axis=-1))

def nearest_neighbors_FV(coords, box, cutoff):
    pbc = True
    x = np.array(coords)
    NNList = []
    for i in range(0,len(coords)):
        NNList.append([])
    for i,pos in enumerate(x):
        dists = distance(pos, x, box, pbc)
        mask = (dists > 1e-15) & (dists <= cutoff)
        for j in mask.nonzero()[0]:
            NNList[i].append(j)
    return NNList

######################################################################
## Generate nearest neighbor list from coords
## Assumes only C present
######################################################################
def nearest_neighbors(coords,abc,cutoff):
    NNList = []
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
                NeighborElements.append(j)
        NNList.append(NeighborElements)
    return(NNList)

######################################################################
## Check if Carbon i is within a 6-membered carbon ring
######################################################################
def carbon_ring(i,NNList,atomtypes):
    for nn in NNList[i]:
        if "C" in atomtypes[nn]:
            for nnn in NNList[nn]:
                if "C" in atomtypes[nnn]:
                    if nnn != i:
                        for nnnn in NNList[nnn]:
                            if "C" in atomtypes[nnnn]:
                                if nnnn != nn:
                                    for nnnnn in NNList[nnnn]:
                                        if "C" in atomtypes[nnnnn]:
                                            if nnnnn != nnn:
                                                for nnnnnn in NNList[nnnnn]:
                                                    if "C" in atomtypes[nnnnnn]:
                                                        if nnnnnn != nnnn:
                                                            for nnnnnnn in NNList[nnnnnn]:
                                                                if nnnnnnn == i:
                                                                    return(True)
    return(False)

######################################################################
## For all atoms, carbon atoms return True/False if they are within 
## a six-membered ring, all other atom types return []
######################################################################
def find_rings(coords, atomtypes, NNList):
    Rings = []
    for i in range(0,len(atomtypes)):
        if "C" in atomtypes[i]:
            RingElement = carbon_ring(i,NNList,atomtypes)
            Rings.append(RingElement)
        else:
            Rings.append([])
    return(Rings)

######################################################################
## Trim carbons not defined by OPLS atom-types
######################################################################
def trim_hanging_carbons(atomtypes, PBC, coords):
    MinDist = 1.5 #C-C bond distance
    newcoords = []
    newatomtypes = []
    numatoms = 0
    moreHangingCarbons = False
    NNList = nearest_neighbors(coords, PBC, MinDist)
    NNNList = nearest_neighbors(coords, PBC, 2.9)
    Rings = find_rings(coords, atomtypes, NNList)
    for i in range(0,len(atomtypes)):
        if len(NNList[i]) == 0:
            pass
        elif Rings[i] == False:
            moreHangingCarbons = True
            pass
        elif ("C" in atomtypes[i] or "O" in atomtypes[i]) and len(NNList[i]) < 2:
            moreHangingCarbons = True
            pass
        #elif ("C" in atomtypes[i] or "O" in atomtypes[i]) and len(NNList[i]) == 2 and ((len(NNNList[(NNList[i][0])]) <= 6) or (len(NNNList[(NNNList[i][1])]) <= 6)):
            #moreHangingCarbons = True
            #pass
        else:
            newcoords.append( coords[i] )
            newatomtypes.append(atomtypes[i])
            numatoms += 1
    return(newcoords, newatomtypes, moreHangingCarbons)

######################################################################
## 
######################################################################
def find_defect_edge_carbons(AtomTypes, PBC, Coords):
    DefectEdgeCarbons = []
    NNNList = nearest_neighbors(Coords, PBC, 2.9)
    for i in range(0,len(AtomTypes)):
        if len(NNNList[i]) < 12:
            DefectEdgeCarbons.append(i)
    return(DefectEdgeCarbons)

######################################################################
## 
######################################################################
def remove_atom(AtomToRemove,AtomTypes,Coords):
    NewAtomTypes = []
    NewCoords = []
    for i in range(0,len(AtomTypes)):
        if i != AtomToRemove:
            NewAtomTypes.append(AtomTypes[i])
            NewCoords.append(Coords[i])
    return(NewAtomTypes,NewCoords)

######################################################################
## TODO update this function
######################################################################
def remove_excess_carbon(AtomTypes, PBC, Coords, OriginalNumCarbons,NumAtomsToRemove):
    print("Removing remaining excess carbon...")
    ## Find all defect edges
    DefectEdgeCarbons = find_defect_edge_carbons(AtomTypes, PBC, Coords)

    ## Remove random defect edge site and re-trim carbons until 
    ## the number of carbons removed = NumAtomsToRemove
    while NumAtomsToRemove > (OriginalNumCarbons - len(Coords)):
        AtomToRemove = random.choice(DefectEdgeCarbons)
        AtomTypes,Coords = remove_atom(AtomToRemove,AtomTypes,Coords)
        moreHangingCarbons = True
        while moreHangingCarbons == True:
            Coords, AtomTypes, moreHangingCarbons = trim_hanging_carbons(AtomTypes, PBC, Coords) 
        
    return(Coords, AtomTypes) 

######################################################################
## Write an output PDB file
######################################################################
def write_pdb_file(coords,atomtypes,abc,FileName):
    newlines = []
    count = 0
    for i in range(0,len(coords)):
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i])))
        count += 1

    file = open(FileName,"w")
    file.write("AUTHOR    GENERATED BY NG TILE TRIMMING\n")
    file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(abc[0],abc[1],abc[2],90.0,90.0,90.0))
    for line in newlines:
        file.write("%s\n"%(line))
    file.close()

    return()

############################
## MAIN			  ##
############################
def main():
    ## Read inputs
    InputFileName,OutFile,NumAtomsToRemove,Orthogonalize = read_inputs()
    ## Read input structure
    Coords,AtomTypes,PBC,OriginAtom = read_pdb_to_structure(InputFileName)
    ## Move origin
    Coords = make_consistant_origin(Coords,OriginAtom)
    ## Convert to square
    if Orthogonalize:
        PBC,Coords = orthogonalize_graphene(PBC,Coords)
        Coords = wrap_orthogonal_coords(Coords,PBC)
    ## Record original number of carbons for output purposes
    OriginalNumCarbons = len(Coords)
    ## Remove overlapping atoms 
    Coords,AtomTypes = remove_overlapping_atoms(Coords,AtomTypes,PBC)
    ## Trim hanging carbons
    i = 0
    print("Trim Cycle %s"%(i))
    Coords, AtomTypes, MoreHangingCarbons = trim_hanging_carbons(AtomTypes, PBC, Coords)
    while MoreHangingCarbons == True:
        i += 1
        print("Trim Cycle %s"%(i))
        Coords, AtomTypes, MoreHangingCarbons = trim_hanging_carbons(AtomTypes, PBC, Coords)
    if NumAtomsToRemove > (OriginalNumCarbons - len(Coords)):
        Coords, AtomTypes = remove_excess_carbon(AtomTypes, PBC, Coords, OriginalNumCarbons,NumAtomsToRemove)

    write_pdb_file(Coords,AtomTypes,PBC,OutFile)
    print("Removed %s carbons of %s"%(OriginalNumCarbons-len(Coords),OriginalNumCarbons))

main()

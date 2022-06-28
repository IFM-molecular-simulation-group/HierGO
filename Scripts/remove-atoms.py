#!/usr/bin/python3
import sys, os, random
import numpy as np
from argparse import ArgumentParser

### A simple code to remove atoms based on given indices and add 
##  hydrogen to new two-fold coordinated carbon atoms
### Usage:
###    python remove-atoms.py --infile InputFileName.pdb --i 1 2 3 --outfile OutputFileName.pdb

######################################################################
## Read pdb to local structure format
######################################################################
def read_pdb_to_structure(InputFile):
   coords = []
   atomtypes = []
   for line in open(InputFile):
       if "CRYST1" in line:
           PBC = line.split()
           abc = [ float(PBC[1]), float(PBC[2]), float(PBC[3]) ]
       elif "ATOM" in line or "HETATM" in line:
           fields = line.split()
           coords.append( [ float(fields[5]), float(fields[6]), float(fields[7]) ] )
           atomtypes.append( fields[2] )

   return(coords,atomtypes,abc)

######################################################################
## Write an output PDB file
######################################################################
def write_pdb_file(coords,atomtypes,abc,FileName):
    newlines = []
    count = 0
    for i in range(0,len(coords)):
        #newlines.append(i)
        #newlines.append(coords[i])
        #print("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i])))
        #newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i])))
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i])))
        count += 1

    print("Writing output file %s"%(FileName))
    file = open(FileName,"w")
    file.write("AUTHOR    GENERATED BY NG TILE DECORATION\n")
    file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(abc[0],abc[1],abc[2],90.0,90.0,90.0))
    for line in newlines:
        file.write("%s\n"%(line))
    file.close()

    return()

######################################################################
######################################################################


######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser(description='A simple code to remove atoms based on given indices and add \
    hydrogen to new two-fold coordinated carbon atoms.')
    parser.add_argument('--infile', default="stiched.pdb", help='input file name (pdb)')
    parser.add_argument('--i', nargs='+', help='<Required> Set flag; indicies of atoms to remove', required=True)
    parser.add_argument('--outfile', nargs='?', default="removed.pdb", type=str, help="output file name (pdb)")
    args = parser.parse_args()
    print("Reading input structure from %s"%(args.infile))
    AtomsToRemove = []
    print("Removing atoms: ",end=" ")
    for Atom in args.i:
        AtomsToRemove.append(int(Atom))
        print(Atom,end=" ")
    print("")
    return(args.infile,AtomsToRemove,args.outfile)

######################################################################
## Remove specified atoms
######################################################################
def remove_atoms(Coords,AtomTypes,AtomsToRemove):
    NewCoords = []
    NewAtomTypes = []
    for i in range(0,len(Coords)):
        AtomToRemove = False
        for b in AtomsToRemove:
            if i == b:
                AtomToRemove = True
        if AtomToRemove == False:
            NewCoords.append(Coords[i])
            NewAtomTypes.append(AtomTypes[i])
    return(NewCoords,NewAtomTypes)

######################################################################
## Generate nearest neighbor list from coords
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
            #if dx >= (abc[0]*0.5):
            #    dx -= abc[0]
            dy = abs(iy - jy)
            #if dy >= (abc[1]*0.5):
            #    dy -= abc[1]
            dz = abs(iz - jz)
            #if dz >= (abc[2]*0.5):
            #    dz -= abc[2]
            d = (dx**2 + dy**2 + dz**2 )**0.5
            if d <= cutoff and i != j:
                NeighborElements.append(j)
        NNList.append(NeighborElements)
    return(NNList)

def calculate_hydrogen_coordinates(c1,c2,c3):
    # move to origin
    ba = np.array(c2) - np.array(c1)
    bc = np.array(c3) - np.array(c1)
    d1 = (ba[0]**2 + ba[1]**2)**0.5
    d2 = (bc[0]**2 + bc[1]**2)**0.5
    ## With SW defects the carbons are not symmetric which causes problems for the 
    ## hydrogen flipping
    if abs(d1 - d2) > 0.1:
        if d1 < d2:
            new = np.average(np.array([ ba, bc ]), axis=0, weights=np.array([2,1]))
        else:
            new = np.average(np.array([ ba, bc ]), axis=0, weights=np.array([1,2]))
    ## Symmetric case is averaged
    else:
        new = np.mean( np.array([ ba, bc ]), axis=0 )
    new = new*np.array([-1,-1])
    theta = abs(np.arctan((new[1]/new[0])))
    new_scaled = np.array([ np.cos(theta), np.sin(theta) ]) #for bond length of 1
    if new[0] < 0:
        new_scaled[0] = new_scaled[0]*-1
    if new[1] < 0:
        new_scaled[1] = new_scaled[1]*-1
    new_scaled += np.array(c1)

    return(new_scaled[0],new_scaled[1])

######################################################################
## Add hydrogens to non-decorated edge carbons
######################################################################
#def add_hydrogen(Coords,ABC,EdgeSites,NNList,NewH):
def add_hydrogen(Coords,atomtypes,ABC,NewH):
    NNList = nearest_neighbors(Coords,ABC,1.75)
    upAtoms = []
    downAtoms = []
    for Atom in range(0,len(Coords)):
        Neighbors = NNList[Atom]
        CNeighbors = []
        for Neighbor in Neighbors:
            if "C" in atomtypes[Neighbor]:
                CNeighbors.append(Neighbor)
        ## carb acid c atoms
        if len(CNeighbors) < 2 and "C" in atomtypes[Atom]:
            #print("< 2 for %s"%(Atom))
            pass
        elif len(Neighbors) < 3 and "C" in atomtypes[Atom]:
            c1x,c1y,c1z = Coords[Atom]
            c2 = [ Coords[CNeighbors[0]][0],Coords[CNeighbors[0]][1] ]
            c3 = [ Coords[CNeighbors[1]][0],Coords[CNeighbors[1]][1] ]
            ## Make sure vectors do not cross boundary
            #c2,c3 = remove_PBC([c1x,c1y],c2,c3,ABC)
            ## Calculate hydrogen coordinates
            hx,hy = calculate_hydrogen_coordinates([c1x,c1y],c2,c3)
            hCoords = [ hx,hy,Coords[Atom][2] ]
            ## Append coordinates
            NewH.append(hCoords)
    return(NewH)

######################################################################
## Add new atoms
######################################################################
def add_new_atoms(coords,atomtypes,NewH):
    for hCoords in NewH:
         #print("H %s %s %s"%(hCoords[0],hCoords[1],hCoords[2]))
         coords.append(hCoords)
         atomtypes.append("H")
    return(coords,atomtypes)

######################################################################
## Add hydrogen to two-fold coordinated carbon
######################################################################
def add_only_hydrogen(coords,atomtypes,abc):
    NewH = add_hydrogen(coords,atomtypes,abc,[])
    coords,atomtypes = add_new_atoms(coords,atomtypes,NewH)
    #newcoords = wrap_coords(coords,abc)
    return(coords,atomtypes)

######################################################################
## MAIN
######################################################################
def main():
    ## Read inputs
    InFile,AtomsToRemove,OutFile = read_inputs()
    ## Read structure
    Coords,AtomTypes,PBC = read_pdb_to_structure(InFile)
    ## Remove specified atoms
    Coords,AtomTypes = remove_atoms(Coords,AtomTypes,AtomsToRemove)
    ## Add hydrogen to unbound atoms
    #Coords,Atomtypes = add_only_hydrogen(Coords,AtomTypes,PBC)
    ## Write file
    write_pdb_file(Coords,AtomTypes,PBC,OutFile)

main()

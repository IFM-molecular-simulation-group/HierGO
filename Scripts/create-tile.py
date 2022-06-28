#!/usr/bin/python3
from argparse import ArgumentParser
import random,sys
import numpy as np
from numpy import array
from numpy.linalg import norm

######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser()

    parser.add_argument('--xfactor', help='size in x dimension, \
    defaults to 20', default=20, type=int)

    parser.add_argument('--yfactor', help='size in y dimension, \
    defaults to 20', default=20, type=int)

    parser.add_argument('--percentvac', default=5, type=float,\
    help='amount of carbon vacancies to add')

    parser.add_argument('--distribution', default='mixed', type=str, \
    choices=['small','medium','large','mixed'])

    ## --keep_hanging_carbon
    parser.add_argument('--keep_hanging_carbon', action='store_true', help='optional \
    flag to keep hanging carbon atoms, which by default are trimmed')

    parser.add_argument('--non_periodic', action='store_true', help='optional \
    flag to remove periodic boundary conditions, which by default are kept')

    ## --outfile filename
    parser.add_argument('--outfile', nargs='?', default="tile.pdb", \
    type=str, help="output file name (pdb)")

    args = parser.parse_args()
    print("Size: %s %s"%(args.xfactor,args.yfactor))
    print("Vacancy Amount: %s"%(args.percentvac))
    print("Distribution: %s"%(args.distribution))
    return(args.outfile,args.xfactor,args.yfactor,args.percentvac,args.distribution,args.keep_hanging_carbon,args.non_periodic)

######################################################################
## Create an orthographic box with a template sheet of C
## Unit cell
##  /
## |
##  \
######################################################################
def create_ortho_template(XFactor,YFactor):
    YRows = int(YFactor*2/3)
    c1 = [ 1.841000,1.775000,0]
    c2 = [0.605000,2.484000,0]
    c3 = [0.605000,3.908000,0]
    c4 = [1.839000,4.620000,0]
    XRowCoords = []
    dx = 2.471
    for i in range(0,XFactor):
        XRowCoords.append([float(c1[0]+dx*i),float(c1[1]),float(c1[2])])
        XRowCoords.append([float(c2[0]+dx*i),float(c2[1]),float(c2[2])])
        XRowCoords.append([float(c3[0]+dx*i),float(c3[1]),float(c3[2])])
        XRowCoords.append([float(c4[0]+dx*i),float(c4[1]),float(c4[2])])
    dy = 4.273000
    Coords = []
    AtomTypes = []

    ## Try and make a square
    if XFactor == YFactor:
        FudgeFactor = int( (((YFactor*4/3)*2.1391) - (XFactor*2.471) ) /dy)
        YRows -= FudgeFactor

    for i in range(0,YRows):
        for b in range(0,len(XRowCoords)):
            Coords.append([float(XRowCoords[b][0]),float(XRowCoords[b][1]+i*dy),float(XRowCoords[b][2])])
            AtomTypes.append("C")
    ABC = [float(XFactor*2.471),float((YRows*2)*2.1391),20.000]
    print("Creating sheet of size %s x %s"%(ABC[0],ABC[1]))
    return(Coords,AtomTypes,ABC)

######################################################################
## Define the minimum and maximum defect sizes based on the type of 
## distribution specified
######################################################################
def define_defect_min_max(DistributionType):
    if DistributionType == "mixed":
        DefectMinMax = [1,30]
    elif DistributionType == "small":
        DefectMinMax = [1,5]
    elif DistributionType == "medium":
        DefectMinMax = [5,10]
    elif DistributionType == "large":
        DefectMinMax = [10,30]
    return(DefectMinMax)

######################################################################
## Choose a random defect size based on the defect min and max
######################################################################
def create_defect_size(DefectMinMax):
    Defect = random.randint(DefectMinMax[0],DefectMinMax[1])
    return(Defect)

######################################################################
## Add 1 to a random element in a list
## (For using up additional vacancies when the number of vacancies is 
##  too low to create an additional vacancy)
######################################################################
def add_1_to_random(DefectList):
    Val = random.choice(range(0,len(DefectList)))
    DefectList[Val] += 1
    return(DefectList)

######################################################################
## Create a list of defects based on the number of atoms to remove and 
## the distribution type. A new defect size is chosen randomly from 
## a size between the minimum and maximum defect size until the 
## remaining number of defect vacancies is less than the maximum 
## allowed defect size, in which case the new defect size is chosen 
## randomly between the minimum defect size and the remaining number 
## of vacancies. If the remaining number of vacancies is less than 
## the minimum defect size, the remaining vacancies are added across 
## the defect list randomly. If this is not possible (presumably 
## because there are no defects in the list) an exception is flagged 
## and the code exits.
######################################################################
def create_defect_sizes(numAtoms,PercentVac,DistributionType):
    AtomsToRemove = int(numAtoms*(PercentVac/100))
    DefectMinMax = define_defect_min_max(DistributionType)
    DefectList = []
    NumDefectAtoms = 0

    while NumDefectAtoms < AtomsToRemove:
        RemainingAtoms = AtomsToRemove - NumDefectAtoms
        ## An Exception for if the remaining atoms cannot fulfill the maximum hole size
        if RemainingAtoms < DefectMinMax[1]:
            ## If the Remaining Atoms fulfill defect conditions, use Remaining Atoms as new max size
            if RemainingAtoms >= DefectMinMax[0]:
                NewDefect = create_defect_size([DefectMinMax[0],RemainingAtoms])
                NumDefectAtoms += NewDefect
                DefectList.append(NewDefect)
            ## If the remaining atoms do not fulfill defect conditions, add the remaining atoms randomly across 
            ## the defect list
            else:
                try:
                    for i in range(1,RemainingAtoms):
                        DefectList = add_1_to_random(DefectList)
                    NumDefectAtoms += RemainingAtoms
                except TypeError:
                    sys.exit("The defect sizes specified are not possible with the template size. Create a larger sheet or use smaller defects.")
        else:
            NewDefect = create_defect_size(DefectMinMax)
            NumDefectAtoms += NewDefect
            DefectList.append(NewDefect)
        #NumDefectAtoms += 1
    print("Creating Defects of Sizes: %s"%(DefectList))
    return(DefectList)

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

def determine_starting_sites(NumAtoms,NNList):
    PossibleDefectSites = []
    for i in range(0,NumAtoms):
        if len(NNList[i]) > 2:
            PossibleDefectSites.append(i)
    return(PossibleDefectSites)

######################################################################
## 
######################################################################
def grow_defect(PossibleDefectSites,DefectSize,Coords,ABC,NNList):
    AtomsToRemove = []
    ## pick first atom
    AtomsToRemove.append( random.choice(PossibleDefectSites) )

    ## subsequent possible sites
    NewPossibleSites = []
    for NN in NNList[AtomsToRemove[0]]:
        NewPossibleSites.append(NN)

    ## loop through the number of atoms to remove
    for i in range(1,DefectSize):
        NewAtom = random.choice(NewPossibleSites)
        AtomsToRemove.append(NewAtom)
        NewPossibleSites.remove(NewAtom)
        for NN in NNList[NewAtom]:
            if NN not in AtomsToRemove:
                NewPossibleSites.append(NN)
    return(AtomsToRemove)
        

######################################################################
## 
######################################################################
def place_defects(Coords,ABC,DefectSizes):
    NNList = nearest_neighbors(Coords,ABC,1.6)
    PossibleDefectSites = determine_starting_sites(len(Coords),NNList)
    AllAtomsToRemove = []
    for DefectSize in DefectSizes:
        AtomsToRemove = grow_defect(PossibleDefectSites,DefectSize,Coords,ABC,NNList)
        for Atom in AtomsToRemove:
            AllAtomsToRemove.append(Atom)
            if Atom in PossibleDefectSites:
                PossibleDefectSites.remove(Atom)
            for NN in NNList[Atom]:
                if NN in PossibleDefectSites:
                    PossibleDefectSites.remove(NN)
                for NNN in NNList[NN]:
                    if NNN in PossibleDefectSites:
                        PossibleDefectSites.remove(NNN)
    
    NewCoords = []
    NewAtomTypes = []
    for i in range(0,len(Coords)):
        if i not in AllAtomsToRemove:
            NewCoords.append(Coords[i])
            NewAtomTypes.append('C')
    return(NewCoords,NewAtomTypes)

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
## Write an output PDB file
######################################################################
def write_pdb_file(coords,atomtypes,abc,FileName):
    newlines = []
    count = 0
    for i in range(0,len(coords)):
        #newlines.append(i)
        #newlines.append(coords[i])
        #print("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2]+31.675,1.00,0.00,str(atomtypes[i])))
        #newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i])))
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i])))
        count += 1

    file = open(FileName,"w")
    file.write("AUTHOR    GENERATED BY NG TILE DECORATION\n")
    file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(abc[0],abc[1],abc[2],90.0,90.0,90.0))
    for line in newlines:
        file.write("%s\n"%(line))
    file.close()

    return()

######################################################################
## MAIN
######################################################################
def main():
    # Read command line arguments
    OutputFileName,XFactor,YFactor,PercentVac,DistributionType,KeepHangingCarbon,NonPeriodic = read_inputs()
    # Returns Template in local format
    Coords,AtomTypes,PBC = create_ortho_template(XFactor,YFactor)
    #print(100 - len(Coords)*(PercentVac/100))
    write_pdb_file(Coords,AtomTypes,PBC,"template.pdb")
    if NonPeriodic:
        PBC[0] += 20
        PBC[1] += 20
        PBC[2] += 20
    # Determine a list of defects which fulfills input conditions
    # (with some exception handling)
    OriginalNumAtoms = len(Coords)
    DefectSizes = create_defect_sizes(OriginalNumAtoms,PercentVac,DistributionType)
    # Place defects
    Coords,AtomTypes = place_defects(Coords,PBC,DefectSizes)
    # Trim 
    if KeepHangingCarbon != True:
        i = 0
        print("Trim Cycle %s"%(i))
        Coords, AtomTypes, MoreHangingCarbons = trim_hanging_carbons(AtomTypes, PBC, Coords)
        while MoreHangingCarbons == True:
            i += 1
            print("Trim Cycle %s"%(i))
            Coords, AtomTypes, MoreHangingCarbons = trim_hanging_carbons(AtomTypes, PBC, Coords)
        ActualPercentVac = 100 - len(Coords)/OriginalNumAtoms*100
        print("Actual vacancy amount of %s percent due to trimming hanging carbon atoms"%(ActualPercentVac))
    else:
        ActualPercentVac = 100 - len(Coords)/OriginalNumAtoms*100
        print("Actual vacancy amount of %s percent"%(ActualPercentVac))
    # Write final tile geometry
    write_pdb_file(Coords,AtomTypes,PBC,OutputFileName)
    #write_output_file(DefectTile,"Tile.cif")

main()

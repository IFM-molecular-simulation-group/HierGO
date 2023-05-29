#!/usr/bin/python3
from argparse import ArgumentParser
import random
import math
import numpy as np

######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser(description='Create a rupture with bridging non-epoxide ether atoms in a graphene lattice.')
    parser.add_argument('--size', nargs='?', default=20, type=int, help='sheet size')
    parser.add_argument('--atoms', nargs='?', default=4, type=int, help='number of C-C atoms to break')
    parser.add_argument('--outfile', nargs='?', default="rupture.pdb", type=str, help='output filename (pdb)')
    args = parser.parse_args()
    return(args.size,args.atoms,args.outfile)

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
## Generate nearest neighbor list from coords
## Assumes only C present
######################################################################
def nearest_neighbor_point(i,coords,abc,cutoff):
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
    return(NeighborElements)

######################################################################
## Create an orthographic box with a template sheet of C
## Unit cell
##  /
## |
##  \
######################################################################
def create_ortho_template(Factor,Square):
    YRows = int(Factor*2/3)
    c1 = [ 1.841000,1.775000,0]
    c2 = [0.605000,2.484000,0]
    c3 = [0.605000,3.908000,0]
    c4 = [1.839000,4.620000,0]
    XRowCoords = []
    dx = 2.471
    for i in range(0,Factor):
        XRowCoords.append([float(c1[0]+dx*i),float(c1[1]),float(c1[2])])
        XRowCoords.append([float(c2[0]+dx*i),float(c2[1]),float(c2[2])])
        XRowCoords.append([float(c3[0]+dx*i),float(c3[1]),float(c3[2])])
        XRowCoords.append([float(c4[0]+dx*i),float(c4[1]),float(c4[2])])
    dy = 4.273000
    Coords = []
    AtomTypes = []

    ## Try and make a square
    if Square:
        FudgeFactor = int( (((Factor*4/3)*2.1391) - (Factor*2.471) ) /dy)
        YRows -= FudgeFactor

    for i in range(0,YRows):
        for b in range(0,len(XRowCoords)):
            Coords.append([float(XRowCoords[b][0]),float(XRowCoords[b][1]+i*dy),float(XRowCoords[b][2])])
            AtomTypes.append("C")
    ABC = [float(Factor*2.471),float((YRows*2)*2.1391),20.000]
    print("Creating sheet of size %s x %s"%(ABC[0],ABC[1]))
    return(Coords,AtomTypes,ABC)
    
######################################################################
## Pick where to locate defect
## Roughly center of the sheet
######################################################################
def pick_defect_site(Coords,ABC):
    x = ABC[0]*0.5
    y = ABC[1]*0.5
    d_min = 1000
    for i in range(0,len(Coords)):
        d = ((x-Coords[i][0])**2 + (y-Coords[i][1])**2)**0.5
        if d < d_min:
            DefectSite = i
            d_min = d
    return(DefectSite)

######################################################################
## Equation of a line in MB form given two carbon atoms
## Goes through the midpoint of the two carbon atoms
######################################################################
def find_perp_line_between_two_atoms(c1,c2):
    ## move points to origin
    dx = (c1[0] + c2[0])*0.5
    dy = (c1[1] + c2[1])*0.5
    c1x = c1[0] - dx
    c1y = c1[1] - dy
    c2x = c2[0] - dx
    c2y = c2[1] - dy
    #print(c1x,c1y,c2x,c2y)
    ## slope of line
    m_inverse = (c2y - c1y)/(c2x-c1x)
    m = -1/m_inverse
    ## move it back
    b = dy - m*dx
    ## for vmd vis
    #print("draw line {%s %0.2f %s} {%0.2f %0.2f %s}"%(0.0,b,0.0,49.4,m*49.4+b,0.0))
    return(m,b)

######################################################################
## 
######################################################################
def find_parallel_atoms(m,b,Coords):
    ParallelAtoms = []
    for i in range(0,len(Coords)):
        ix,iy = Coords[i][0],Coords[i][1]
        
        if (abs(m) < 0.01) and (abs((m*ix + b)-iy) < 0.8):
            if i not in ParallelAtoms:
                ParallelAtoms.append(i)
        elif (abs(m) > 0.01) and (abs((m*ix + b)-iy) < 1.7):
            if i not in ParallelAtoms:
                ParallelAtoms.append(i)
    #print("ParallelAtoms: ")
    #for i in ParallelAtoms:
    #    print(i,end=" ")
    #print("")
    return(ParallelAtoms)

######################################################################
## 
######################################################################
def group_parallel_atoms(ParallelAtoms,coords,abc):
    ParallelAtomsGrouped = []
    ParallelAtomsRemaining = []

    for i in ParallelAtoms:
        ParallelAtomsRemaining.append(i)

    for i in ParallelAtoms:
        if i in ParallelAtomsRemaining:
            ## could save time by only looping over parallel atom coords
            NNList = nearest_neighbor_point(i,coords,abc,1.8)
            for nn in NNList:
                if nn in ParallelAtomsRemaining:
                    ParallelAtomsGrouped.append([i,nn])
                    ParallelAtomsRemaining.remove(i)
                    ParallelAtomsRemaining.remove(nn)
    
    return(ParallelAtomsGrouped)

######################################################################
## 
######################################################################
def find_bonds_parallel(Atoms,Coords,ABC,NNList):
    #print("Initial Atoms:")
    #print(Atoms[0],Atoms[1])
    m,b = find_perp_line_between_two_atoms(Coords[Atoms[0]],Coords[Atoms[1]])
    ParallelAtoms = find_parallel_atoms(m,b,Coords)
    ParallelAtomsGrouped = group_parallel_atoms(ParallelAtoms,Coords,ABC)
    return(ParallelAtomsGrouped)

######################################################################
## 
######################################################################
def choose_defect_carbons(BondList,ParallelAtomsGrouped,Coords,RuptureSize):
    BondListAtoms = []
    for i in BondList:
        BondListAtoms.append(i[0])
        BondListAtoms.append(i[1])
    for i in range(0,RuptureSize):
        PossibleNewBond = []
        for ParallelAtoms in ParallelAtomsGrouped:
            if ParallelAtoms[0] not in BondListAtoms:
                for OldBond in BondList:
                    dx = abs(Coords[ParallelAtoms[0]][0] - Coords[OldBond[0]][0])
                    dy = abs(Coords[ParallelAtoms[0]][1] - Coords[OldBond[0]][1])
                    d1 = (dx**2 + dy**2)**0.5
                    dx1 = abs(Coords[ParallelAtoms[0]][0] - Coords[OldBond[1]][0])
                    dy1 = abs(Coords[ParallelAtoms[0]][1] - Coords[OldBond[1]][1])
                    d2 = (dx1**2 + dy1**2)**0.5
                    if d1 < 2.5 or d2 < 2.5:
                        PossibleNewBond.append(ParallelAtoms)
        NewBond = random.choice(PossibleNewBond)
        BondList.append(NewBond)
        BondListAtoms.append(NewBond[0])
        BondListAtoms.append(NewBond[1])
    return(BondList)

def add_oxygen(c1,c2,up):
    ox = (c1[0]+c2[0])*0.5
    oy = (c1[1]+c2[1])*0.5
    z = (c1[2]+c2[2])*0.5
    if up == True:
        oz = z + 0.6
    else:
        oz = z - 0.6
    oCoords = [ox,oy,oz]
    return(oCoords)

######################################################################
## 
######################################################################
def break_bonds_add_oxygen(BondList,Coords,AtomTypes):
    ## reorder bondlist by avg x value for oxygen z direction
    OrderedBondList = []
    while len(OrderedBondList) < len(BondList):
        dx_min = 100000
        for Bond in BondList:
            if Bond not in OrderedBondList:
                c1x = Coords[Bond[0]][0]
                c2x = Coords[Bond[1]][0]
                dx = (c2x + c1x)*0.5
                if dx <= dx_min:
                    dx_min = dx
                    NewBond = Bond
        OrderedBondList.append(NewBond)

    up = True
    for Bond in OrderedBondList:
        c1x = Coords[Bond[0]][0]
        c1y = Coords[Bond[0]][1]
        c2x = Coords[Bond[1]][0]
        c2y = Coords[Bond[1]][1]
        dx = (c2x - c1x)*0.3
        dy = (c2y - c1y)*0.3
        if c1x < c2x:
            c1xnew = c1x - dx
            c2xnew = c2x + dx
        else:
            c1xnew = c1x + dx
            c2xnew = c2x - dx
        if c1y < c2y:
            c1ynew = c1y - dy
            c2ynew = c2y + dy
        else:
            c1ynew = c1y + dy
            c2ynew = c2y - dy
        if (c1x > c2x and c1y < c2y):
            c1xnew = c1x + abs(dx)
            c2xnew = c2x - abs(dx)
            c1ynew = c1y - abs(dy)
            c2ynew = c2y + abs(dy)
        elif (c1x < c2x and c1y > c2y):
            c1xnew = c1x - abs(dx)
            c2xnew = c2x + abs(dx)
            c1ynew = c1y + abs(dy)
            c2ynew = c2y - abs(dy)
        Coords[Bond[0]] = [ c1xnew, c1ynew, Coords[Bond[0]][2] ]
        Coords[Bond[1]] = [ c2xnew, c2ynew, Coords[Bond[1]][2] ]
        if up == True:
            oCoords = add_oxygen(Coords[Bond[0]],Coords[Bond[1]],True)
            up = False
        else:
            oCoords = add_oxygen(Coords[Bond[0]],Coords[Bond[1]],False)
            up = True
        Coords.append(oCoords)
        AtomTypes.append("O")
    return(Coords,AtomTypes)

######################################################################
## Centered at the defect site, break C-C bonds by elongation and 
## add ether atoms across the rupture with alternating z
######################################################################
def create_rupture(Coords,ABC,AtomTypes,DefectSite,RuptureSize):
    print("Creating rupture of size %s"%(RuptureSize))
    ## first bond to lengthen is randomly chosen from NN
    NNList = nearest_neighbors(Coords,ABC,1.7)
    BondList = [ [DefectSite,random.choice(NNList[DefectSite])] ]
    if Coords[BondList[0][1]][1] < Coords[BondList[0][0]][1] and Coords[BondList[0][1]][0] < Coords[BondList[0][0]][0]:
        BondList = [ [BondList[0][1],BondList[0][0]] ]
    ## for now use chosen sites
    #BondList = [ [359,436] ]

    ## find all parallel atoms
    ParallelAtomsGrouped = find_bonds_parallel(BondList[0],Coords,ABC,NNList)

    ## Choose atoms
    BondList = choose_defect_carbons(BondList,ParallelAtomsGrouped,Coords,(RuptureSize-1))
    #print("DefectAtoms: ")
    #for i in BondList:
    #    print(i[0],i[1],end=" ")
    #print("")

    ## Lengthen C-C bond and add oxygen
    Coords,AtomTypes = break_bonds_add_oxygen(BondList,Coords,AtomTypes)

    return(Coords,AtomTypes)

######################################################################
## Write an output PDB file
######################################################################
def write_pdb_file(Coords,AtomTypes,ABC,FileName):
    newlines = []
    count = 0
    for i in range(0,len(Coords)):
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,AtomTypes[i],"UNL",1,Coords[i][0],Coords[i][1],Coords[i][2],1.00,0.00,AtomTypes[i]))
        count += 1

    print("Writing output file %s"%(FileName))
    file = open(FileName,"w")
    file.write("AUTHOR    GENERATED BY NG TILE DECORATION\n")
    #file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(InputStructure.lattice.abc[1],InputStructure.lattice.abc[2],InputStructure.lattice.abc[0],InputStructure.lattice.beta,InputStructure.lattice.gamma,InputStructure.lattice.alpha))
    file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(ABC[0],ABC[1],ABC[2],90,90,90))
    for line in newlines:
        file.write("%s\n"%(line))
    file.close()

    return()

######################################################################
## MAIN
######################################################################
def main():
    Square = True
    Factor,RuptureSize,OutputFile =  read_inputs()
    Coords,AtomTypes,ABC = create_ortho_template(Factor,Square)
    #write_pdb_file(Coords,AtomTypes,ABC,"template.pdb")
    DefectSite = pick_defect_site(Coords,ABC)
    Coords,AtomTypes = create_rupture(Coords,ABC,AtomTypes,DefectSite,RuptureSize)
    write_pdb_file(Coords,AtomTypes,ABC,OutputFile) 

main()

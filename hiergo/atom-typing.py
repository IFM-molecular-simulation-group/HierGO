#!/usr/bin/python3
from argparse import ArgumentParser
import random,sys
import numpy as np
import time
from collections import Counter

#-infile ortho.pdb -outfile typed.pdb
## TODO: flag multiple functional groups on the same C atom 

######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser(description='Type atoms according to OPLS atom types.')
    parser.add_argument('--infile', help='<Required> Set flag; input file name (pdb)', required=True)
    parser.add_argument('--outfile', nargs='?', default="typed.pdb", type=str, help="output file name (pdb)")
    parser.add_argument('--pbc', action='store_true', help='default to ignoring PBC')

    args = parser.parse_args()
    print("Reading input structure from %s"%(args.infile))
    return(args.infile,args.outfile,args.pbc)

######################################################################
## Read pdb to local structure format
######################################################################
def read_pdb_to_structure(InputFile):
   Coords = np.array([])
   AtomTypes = []
   for line in open(InputFile):
       if "CRYST1" in line:
           PBC = line.split()
           ABC = np.array([ float(PBC[1]), float(PBC[2]), float(PBC[3]) ])
       elif "ATOM" in line or "HETATM" in line:
           fields = line.split()
           #coords.append( [ float(fields[5]), float(fields[6]), float(fields[7]) ] )
           if Coords.size == 0:
               Coords = np.append(Coords, [ float(fields[5]), float(fields[6]), float(fields[7]) ])
           else:
               Coords = np.vstack([Coords, [ float(fields[5]), float(fields[6]), float(fields[7]) ]])
           AtomTypes.append( fields[2] )

   return(Coords,AtomTypes,ABC)

######################################################################
## Generate nearest neighbor list from coords
######################################################################
def find_nearest_neighbors(coords,atomtypes,abc):
    CNNList = []
    ONNList = []
    HNNList = []
    for i in range(0,len(coords)):
        #print("Atom %s"%(i))
        if atomtypes[i] == "C":
            CCutoff = 1.7
            OCutoff = 1.81
            HCutoff = 1.3
        elif atomtypes[i] == "O":
            CCutoff = 1.81
            OCutoff = 1.81
            HCutoff = 1.3
        elif atomtypes[i] == "H":
            CCutoff = 1.3
            OCutoff = 1.3
            HCutoff = 1.3
        ix = coords[i][0]
        iy = coords[i][1]
        iz = coords[i][2]
        CNeighborElements = []
        ONeighborElements = []
        HNeighborElements = []
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
            if d <= CCutoff and i != j and atomtypes[j] == "C":
                CNeighborElements.append(j)
            elif d <= OCutoff and i != j and atomtypes[j] == "O":
                ONeighborElements.append(j)
            elif d <= HCutoff and i != j and atomtypes[j] == "H":
                HNeighborElements.append(j)
        CNNList.append(CNeighborElements)
        ONNList.append(ONeighborElements)
        HNNList.append(HNeighborElements)
    return(CNNList,ONNList,HNNList)

######################################################################
## FV NN
######################################################################
def distance(x0, x1, box, pbc):
    # xo is a position of one atom, x1 is an array of positions
    # use the pbc bool mask to set the periodicity
    delta = np.abs(x0 - x1)
    delta[:,pbc] -= box[pbc] * np.round(delta[:,pbc]/(box[pbc]))
    return(np.sqrt((delta ** 2).sum(axis=-1)))

def nearest_neighbors_fv(Coords, box, Cutoff):
    pbc = np.full(3, True)
    NNList = []
    for i in range(0,len(Coords)):
        NNList.append([])
    for i,pos in enumerate(Coords):
        dists = distance(pos, Coords, box, pbc)
        mask = (dists > 1e-15) & (dists <= Cutoff)
        for j in mask.nonzero()[0]:
            NNList[i].append(j)
    return(NNList)

def neighbors_by_atomtype(coords,atomtypes,abc):
    NNList = nearest_neighbors_fv(coords, abc, 1.81)
    CNNList = []
    ONNList = []
    HNNList = []
    for i in range(0,len(coords)):
        if atomtypes[i] == "C":
            CCutoff = 1.7
            OCutoff = 1.81
            HCutoff = 1.3
        elif atomtypes[i] == "O":
            CCutoff = 1.81
            OCutoff = 1.81
            HCutoff = 1.3
        elif atomtypes[i] == "H":
            CCutoff = 1.3
            OCutoff = 1.3
            HCutoff = 1.3
        ix = coords[i][0]
        iy = coords[i][1]
        iz = coords[i][2]
        CNeighborElements = []
        ONeighborElements = []
        HNeighborElements = []
        for j in NNList[i]:
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
            if d <= CCutoff and i != j and atomtypes[j] == "C":
                CNeighborElements.append(j)
            elif d <= OCutoff and i != j and atomtypes[j] == "O":
                ONeighborElements.append(j)
            elif d <= HCutoff and i != j and atomtypes[j] == "H":
                HNeighborElements.append(j)
        CNNList.append(CNeighborElements)
        ONNList.append(ONeighborElements)
        HNNList.append(HNeighborElements)
    return(CNNList,ONNList,HNNList)

######################################################################
#opls_155   HYD  HH   0.4180  ###hydroxyl H
#????????   HYD OH1   ??????  ###hydroxyl O
#opls_159   HYD  CD   0.2650  ###hydroxyl C
#opls_180   EPX  OE  -0.4000  ###epoxy    O
#opls_184   EPX  CE   0.2000  ###epoxy    C  
#opls_180   ETH  OT  -0.2300  ###ether    O
#opls_518   ETH  CT   0.0850  ###ether    C
#opls_281   CAB  ON  -0.4700  ###carbonyl O
#opls_280   CAB  CN   0.4700  ###carbonyl C
#opls_269  CARB  OC  -0.4400  ###carboxylic OC
#opls_268  CARB  OH  -0.5300  ###carboxylic OH
#opls_270  CARB HC1   0.4500  ###carboxylic H
#opls_267  CARB  CA   0.5200  ###carboxylic C
#opls_147  CARB  CC   0.0000  ###C-carboxylic
#opls_146   GRA  HC   0.1150  ###terminal H
#opls_145   GRA  CH  -0.1150  ###terminal C
#opls_147   GRA   C   0.0000  ###In-plane sp2 carbon
######################################################################
def opls_atom_typing(coords,atomtypes,abc):
    CNNList,ONNList,HNNList = neighbors_by_atomtype(coords,atomtypes,abc)
    #CNNList,ONNList,HNNList = find_nearest_neighbors(coords,atomtypes,abc)
    newatomtypes = []
    for atom in atomtypes:
        newatomtypes.append(atom)
    for i in range(0,len(coords)):
        if atomtypes[i] == "C":
            if len(CNNList[i]) == 3 and len(ONNList[i]) == 0 and len(HNNList[i]) == 0:
                if len(CNNList[CNNList[i][0]]) == 1 or len(CNNList[CNNList[i][1]]) == 1 or len(CNNList[CNNList[i][2]]) == 1:
                    newatomtypes[i] = "CC"
                else:
                    newatomtypes[i] = "C"
            elif len(CNNList[i]) == 1: #the only hanging carbons should be carboxylic acid groups
                newatomtypes[i] = "CA"
            elif len(ONNList[i]) == 0 and len(HNNList[i]) > 0:
                newatomtypes[i] = "CH"
            elif len(ONNList[i]) == 1 and len(HNNList[ONNList[i][0]]) == 0:
                #print(CNNList[ONNList[i][0]])
                #print("ok")
                #print(CNNList[CNNList[ONNList[i][0]][1]])
                if len(CNNList[ONNList[i][0]]) == 1:
                    ## TODO temp fix of carbonyls being added 
                    newatomtypes[i] = "CN"
                    #newatomtypes[i] = "C"
                elif CNNList[ONNList[i][0]][0] in CNNList[CNNList[ONNList[i][0]][1]]:
                    newatomtypes[i] = "CE"
                else:
                    newatomtypes[i] = "CT"
            elif len(ONNList[i]) == 1 and len(HNNList[ONNList[i][0]]) == 1:
                newatomtypes[i] = "CD"
        elif atomtypes[i] == "O":
            if len(CNNList[i]) == 2 and len(HNNList[i]) == 0:
                if CNNList[i][0] in CNNList[CNNList[i][1]]:
                    newatomtypes[i] = "OE"
                else:
                    newatomtypes[i] = "OT"
            elif len(CNNList[i]) == 1 and len(HNNList[i]) == 1:
                if len(ONNList[CNNList[i][0]]) == 2:
                    newatomtypes[i] = "OH"
                    ## hydroxyls inadvertedly bonded to epoxy groups can be mislabelled as carboxylic OH groups, remove these
                    if (len(CNNList[ONNList[CNNList[i][0]][0]]) > 1) or (len(CNNList[ONNList[CNNList[i][0]][1]]) > 1):
                        newatomtypes[i] = "OU"
                else:
                    newatomtypes[i] = "OH1"
            elif len(CNNList[i]) == 1 and len(HNNList[i]) == 0: 
                if len(ONNList[CNNList[i][0]]) == 2:
                    newatomtypes[i] = "OC"
                else:
                    ## TODO: temp fix of carbonyls being added
                    newatomtypes[i] = "ON"
                    #newatomtypes[i] = "OU"
            elif len(CNNList[i]) == 0:
                newatomtypes[i] = "OU" #unbound oxygen, will get removed
            else:
                print("An oxygen atom was unable to be typed. Removing atom %s. Is geometry sensible? Overlapping atoms will cause problems."%(i))
        elif atomtypes[i] == "H":
            if len(CNNList[i]) > 0:
                newatomtypes[i] = "HC"
            elif len(ONNList[i]) > 0:
                newatomtypes[i] = "HU" #unbound hydrogen, will get removed
                for o in ONNList[i]:
                    if len(CNNList[o]) > 0:
                        if len(ONNList[CNNList[o][0]]) == 2:
                            #print("Atom: %s, ONNList[i]: %s, CNNList[o]: %s, ONNList[CNNList[o][0]]: %s"%(i,ONNList[i],CNNList[o],ONNList[CNNList[o][0]]))
                            ##TODO: fix this, will work for epoxy/hydrox
                            #newatomtypes[i] = "HU"
                            newatomtypes[i] = "HC1"
                        else:
                            newatomtypes[i] = "HH"
            else:
                newatomtypes[i] = "HU"
                print("A hydrogen atom was unable to be typed. Removing atom %s. Is geometry sensible? Overlapping atoms will cause problems."%(i))
    return(newatomtypes)

######################################################################
# Remove unbound oxygen and hydrogen atoms
# (assumes no unbound carbon)
######################################################################
def remove_unbound_atoms(coords,atomtypes):
    newcoords = np.array([])
    newatomtypes = []
    UnboundAtoms = []
    for i in range(0,len(coords)):
        if atomtypes[i] != "O" and atomtypes[i] != "H" and atomtypes[i] != "OU" and atomtypes[i] != "HU":
            if newcoords.size == 0:
                newcoords = np.append(newcoords,[coords[i]])
            else:
                newcoords = np.vstack([newcoords,coords[i]])
            newatomtypes.append(atomtypes[i])
        else:
            UnboundAtoms.append(i)
    if UnboundAtoms != []:
        print("Unbound atoms removed: ")
        for i in UnboundAtoms:
            print(i,end=" ")
        print("")
    return(newcoords,newatomtypes)

######################################################################
## Generate nearest neighbor list from coords
######################################################################
def nearest_neighbors(coords,abc,atomtypes):
    NNList = []
    for i in range(0,len(coords)):
        if atomtypes[i][:1] == "C":
            CCutoff = 1.7
            OCutoff = 1.81
            HCutoff = 1.3
        elif atomtypes[i][:1] == "O":
            CCutoff = 1.7
            OCutoff = 1.81
            HCutoff = 1.3
        elif atomtypes[i][:1] == "H":
            CCutoff = 1.3
            OCutoff = 1.3
            HCutoff = 1.3
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
            if d <= CCutoff and i != j and atomtypes[j][:1] == "C":
                NeighborElements.append(j)
            elif d <= OCutoff and i != j and atomtypes[j][:1] == "O":
                NeighborElements.append(j)
            elif d <= HCutoff and i != j and atomtypes[j][:1] == "H":
                NeighborElements.append(j)
        Neighbor_String = ("CONECT{:>5s}".format(str(i)))
        #print(NeighborElements)
        #print(NeighborElements.sorted())
        for Neighbor in sorted(NeighborElements):
            #{:6s}{:5d} {:^4s} {:3s}  {:4d}
            Neighbor_String += ("{:>5s}".format(str(Neighbor)))
        Neighbor_String += "\n"
        NNList.append(Neighbor_String)
    return(NNList)

######################################################################
#opls_155   HYD  HH   0.4180  ###hydroxyl H
#????????   HYD OH1   ??????  ###hydroxyl O
#opls_159   HYD  CD   0.2650  ###hydroxyl C
#opls_180   EPX  OE  -0.4000  ###epoxy    O
#opls_184   EPX  CE   0.2000  ###epoxy    C  
#opls_180   ETH  OT  -0.2300  ###ether    O
#opls_518   ETH  CT   0.0850  ###ether    C
#opls_281   CAB  ON  -0.4700  ###carbonyl O
#opls_280   CAB  CN   0.4700  ###carbonyl C
#opls_269  CARB  OC  -0.4400  ###carboxylic OC
#opls_268  CARB  OH  -0.5300  ###carboxylic OH
#opls_270  CARB HC1   0.4500  ###carboxylic H
#opls_267  CARB  CA   0.5200  ###carboxylic C
#opls_147  CARB  CC   0.0000  ###C-carboxylic
#opls_146   GRA  HC   0.1150  ###terminal H
#opls_145   GRA  CH  -0.1150  ###terminal C
#opls_147   GRA   C   0.0000  ###In-plane sp2 carbon
######################################################################
def check_typing_numbers(atomtypes):
    types = Counter(atomtypes)
    warnings = 0
    if types['CE']/2 != types['OE']:
        print("Warning: %s carbon epoxy atoms when there should be %s"%(types['CE'],types['OE']*2))
        warnings += 1
    if types['OH1'] != types['HH'] or types['OH1'] != types['CD']:
        print("Warning: %s hydroxyl oxygen, %s hydroxyl hydrogen and %s hydroxyl carbon"%(types['OH1'],types['HH'],types['CD']))
        warnings += 1
    if types['CT']/2 != types['OT']:
        print("Warning: %s carbon ether atoms when there should be %s"%(types['CT'],types['OT']*2))
        warnings += 1
    if types['CN'] != types['ON']:
        print("Warning: %s carbonyl carbon atoms and %s carbonyl oxygen atoms"%(types['CN'],types['ON']))
        warnings += 1
    if types['CA'] != types['OH'] or types['CA'] != types['OC'] or types['CA'] != types['HC1'] or types['CA'] != types['CC']:
        print("Warning: %s CA, %s OH, %s OC, %s HC1, %s CC"%(types['CA'],types['OH'],types['OC'],types['HC1'],types['CC']))
        warnings += 1
    if types['HC'] != types['CH']:
        print("Warning: %s HC, %s CH"%(types['HC'],types['CH']))
        warnings += 1

    if warnings == 0:
        print("Passed typing equivalency check.")
    return()

######################################################################
# print xyz coords
######################################################################
def write_xyz_file(coords,atomtypes,pbc,FileName):
    file = open(FileName,"w")
    file.write("%s\n"%(len(coords)))
    file.write("%s %s %s\n"%(pbc[0],pbc[1],pbc[2]))
    for i in range(0,len(coords)):
        file.write("%s %s %s %s\n"%(atomtypes[i],coords[i][0],coords[i][1],coords[i][2]))
    file.close()
    return()

######################################################################
## Write an output PDB file
######################################################################
def write_pdb_file(coords,atomtypes,abc,FileName):
    newlines = []
    count = 0
    for i in range(0,len(coords)):
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n".format("ATOM",count,str(atomtypes[i]),"UNL",1,coords[i][0],coords[i][1],coords[i][2],1.00,0.00,str(atomtypes[i][:1])))
        count += 1

    ## CONECT List 
    #for i in range(0,len(coords)):
    #    if NNList[i] != "":
    #        newlines.append(NNList[i])

    print("Writing output file %s"%(FileName))
    file = open(FileName,"w")
    file.write("AUTHOR    GENERATED BY NG TILE TYPING\n")
    file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(abc[0],abc[1],abc[2],90.0,90.0,90.0))
    for line in newlines:
        file.write("%s"%(line))
    file.close()

    return()

######################################################################
## MAIN
######################################################################
def main():
    t0 = time.perf_counter()
    InputFile,OutputFile,Periodic = read_inputs()
    coords,atomtypes,abc = read_pdb_to_structure(InputFile)
    if Periodic != True:
       abc[0] += 20
       abc[1] += 20
       abc[2] += 20
    #write_pdb_file(coords,atomtypes,abc,"ordered.pdb")
    #sys.exit()
    #coords,atomtypes = remove_overlapping_atoms(coords,atomtypes,abc)
    atomtypes = opls_atom_typing(coords,atomtypes,abc)
    coords,atomtypes = remove_unbound_atoms(coords,atomtypes)
    #NNList = nearest_neighbors_fv(coords,abc,atomtypes)
    #write_xyz_file(coords,atomtypes,abc,"typed.xyz")
    write_pdb_file(coords,atomtypes,abc,OutputFile)
    check_typing_numbers(atomtypes)
    t1 = time.perf_counter()
    print(f'Result computed in {t1 - t0:0.4f} seconds')
main()



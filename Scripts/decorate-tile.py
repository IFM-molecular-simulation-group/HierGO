#!/usr/bin/python3
from argparse import ArgumentParser
import random,sys
import numpy as np
from numpy import array
from numpy.linalg import norm
from collections import Counter

#-infile ../TileCreation/20-20-mixed-oredo.pdb -percento 10 -otype 1:1
#python ../Scripts/decorate-tile-2.py -infile 20-10-mixed-redo-o-sw.pdb -percento 25 -otype 1:1

## LIMITATIONS
#  - assumes flat input structure

######################################################################
## General Mathmatical functions
######################################################################

def get_unit_normal(a, b):
    """
    Takes two vectors and returns a unit normal vector with vector 'a' as the
    origin
    """

    # centre vector about <0,0>
    a_ = array([0, 0])
    b_ = b - a

    norm_vector = array([-b_[1], b_[0]])
    norm_vector_o = array([b_[1], -b_[0]])

    unit_norm_vector = norm_vector / norm(
        norm_vector)  # <-b_[1], b_[0]> / |<-b_[1], b_[0]>|

    unit_norm_vector_o = norm_vector_o / norm(
        norm_vector_o)  # <-b_[1], b_[0]> / |<-b_[1], b_[0]>|

    return unit_norm_vector, unit_norm_vector_o


def scale_unit_vector(unit_vector, length):
    """
    scales the unit vector by the length parameter
    """
    return length * unit_vector


def midpoint_vector(vect_a, vect_b):
    """
    returns a vector to the mid point of a segment
    """
    from numpy import array

    return array([(vect_b[0] + vect_a[0]) / 2
                     , (vect_b[1] + vect_a[1]) / 2
                  ])


def locate_vetor(in_vector, location_vector):
    """
    locates a vector with origin <0,0> at some new location_vector
    """
    return in_vector + location_vector

def get_new_position(a, b, length):
    """
    Takes two vectors and returns a new postion normal to the segment between a
    and b, with given length and located at the midpoint.

    eg input.
    a = array([0, 1])
    b = array([2, 3])

    :param a: np.array vector type for postion a
    :param b: np.array vector type for postion b
    :param length: float for scaling length

    :returns: np.array vector type for new position
    """
    # Calc. unit normal vector about a
    unit_normal_vect,unit_normal_vect_o = get_unit_normal(a, b)

    # scale the unit vector to the right length
    scaled_unit_vect = scale_unit_vector(
        unit_vector=unit_normal_vect
        , length=length
    )

    scaled_unit_vect_o = scale_unit_vector(
        unit_vector=unit_normal_vect_o
        , length=length
    )

    # Calc. the midpoint vector
    midpoint_vect = midpoint_vector(a, b)

    # locate vector at the midpoint
    position_vect = locate_vetor(
        in_vector=scaled_unit_vect
        , location_vector=midpoint_vect
    )
    position_vect_o = locate_vetor(
        in_vector=scaled_unit_vect_o
        , location_vector=midpoint_vect
    )

    return {'a': a
        , 'b': b
        , 'length': length
        , 'unit_normal_vect': unit_normal_vect
        , 'scaled_unit_vect': scaled_unit_vect
        , 'midpoint_vect': midpoint_vect
        , 'position_vect': position_vect
        , 'position_vect_o': position_vect_o
            }

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy

######################################################################
## Decide top or bottom randomly
######################################################################
def coin_flip():
    Side = random.choice([True,False])
    if Side == True:
        return(True)
    return(False)

######################################################################
## Read inputs
## TODO: add aged/fresh and for fresh epoxy:hydroxyl ratio, whereas 
##       for aged a temperature
##       use edges first y/n
######################################################################
def read_inputs():
    parser = ArgumentParser(description='A script which adds oxygen content \
    of a given proportion based on input arguments. Ether the --fresh or the \
    --aged flag must be specified.')

    ## --infile filename
    parser.add_argument('--infile', help='input file name (pdb)', \
    default="tile.pdb")

    ## --outfile filename
    parser.add_argument('--outfile', nargs='?', default="dec-tile.pdb", \
    type=str, help="output file name (pdb)")

    ## --percento 30
    parser.add_argument('--percento', help='overall amount of oxygen, \
    oxygen will be added to the edges in order to satisfy this percent, \
    defaults to 0 which will just add hydrogen', default=0, type=float)

    group = parser.add_mutually_exclusive_group()
    ## --fresh 1:1
    group.add_argument('--fresh', nargs='?', const="1:1", type=str, \
    help="decorate with hydroxyl and epoxy groups only in a ratio \
    specified by hydroxyl:epoxy, defaults to 1:1.")

    ## --aged 900
    group.add_argument('--aged', nargs='?', const=900, type=int, \
    choices=[600,900,1200,2500], help='decorate with hydroxyl, epoxy, \
    carboxylic acid, epoxide and non-epoxide ether in a ratio specified \
    by an annealing ``temperature.\'\' Defaults to 900K.')

    ## --bunching 5 10
    ## --bunching 15,15,15 10
    parser.add_argument('--bunching', nargs='+', type=str, help='first place oxygen \
    in a circle of a given diameter around a given atom index or coordinate. \
    Syntax --bunching [Atom #] [Diameter]')

    ## --decorate-edges-first
    parser.add_argument('--decorate_edges_first', action='store_true', help='optional \
    flag with no arguments to decorate defect edges before the basal plane')

    ## --exclude_tile_edges
    parser.add_argument('--exclude_tile_edges', action='store_true', help='optional \
    flag with no arguments to exclude tile edges in the decoration process')

    ## print statements and further formatting of parsed arguments
    args = parser.parse_args()
    profile = None
    bunch_origin = None
    bunch_diameter = None
    fresh = None

    print("Input parameters:")
    print("    Reading file %s"%(args.infile))
    print("    %s%% oxygen"%(args.percento))

    if args.fresh != None:
        elements = args.fresh.split(":")
        profile = [ int(elements[0]), int(elements[1]) ]
        print("    [Epoxy,Hydroxyl] = %s"%(args.fresh))
        fresh = True
    elif args.aged != None:
        print("    Aging profile of %sK"%(args.aged))
        profile = args.aged
        fresh = False

    if args.bunching != None:
        bunch_diameter = float(args.bunching[1])
        if "," in args.bunching[0]:
            bunch_origin = args.bunching[0].split(",")
            print("    Placing oxygen in a circle of diameter %s centered around %s,%s,%s"%(bunch_diameter,bunch_origin[0],bunch_origin[1],bunch_origin[2]))
        else:
            bunch_origin = int(args.bunching[0])
            print("    Placing oxygen in a circle of diameter %s centered around atom %s"%(bunch_diameter,bunch_origin))

    return(args.infile,args.outfile,args.percento,fresh,profile,bunch_origin,bunch_diameter,args.decorate_edges_first,args.exclude_tile_edges)

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



def determine_oxygen(atomtypes):
    c = Counter(atomtypes)
    OriginalO = 0
    OriginalC = 0
    OriginalH = 0
    for i in c:
        if "C" in i:
            OriginalC = c[i]
        elif "O" in i:
            OriginalO = c[i]
        elif "H" in i:
            OriginalH = c[i]
        
    return(OriginalO,OriginalC,OriginalH)

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
######################################################################
def nearest_neighbor_single(i,coords,abc,cutoff):
    ix,iy,iz = i[0],i[1],i[2]
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
        if d <= cutoff and abs(d) > 0.001:
            NeighborElements.append(j)
    return(NeighborElements)

######################################################################
## Calculate number of total oxygens, num epoxys and num hydroxyls 
## and return in dict
######################################################################
def calculate_oxygen_num_type(NumC,ONumO,PercentO,HERatio):
    if ONumO != 0:
        print("Reading structure with %0.2f%% initial oxygen"%(ONumO/NumC*100))
    NumO = int(NumC * (PercentO / 100 )) - ONumO
    HFrac = HERatio[0] / (HERatio[0] + HERatio[1])
    EFrac = HERatio[1] / (HERatio[0] + HERatio[1])
    NumHy = int(HFrac*NumO)
    NumEp = int(EFrac*NumO)

    OxygenNumAndType = {}
    # 30% hydroxyl; 70% dihydroxyls
    OxygenNumAndType["Surface_Hydroxyl"] = int(round(NumHy*0.3,0))
    OxygenNumAndType["Dihydroxyl"] = int(round(NumHy*0.7,0))
    if (OxygenNumAndType["Dihydroxyl"]%2) != 0:
        OxygenNumAndType["Dihydroxyl"] -= 1
        OxygenNumAndType["Surface_Hydroxyl"] += 1
    OxygenNumAndType["Epoxy"] = NumEp
    OxygenNumAndType["Total"] = NumO

    print("Adding %s oxygen atoms:"%(OxygenNumAndType["Total"]))
    print("       O_hydroxyl:   %s"%(OxygenNumAndType["Surface_Hydroxyl"]))
    print("     O_dihydroxyl:   %s"%(OxygenNumAndType["Dihydroxyl"]))
    print("          O_epoxy:   %s"%(OxygenNumAndType["Epoxy"]))
    #print("    O_epoxy"%(NumO,NumHy,NumEp))

    return(OxygenNumAndType)

######################################################################
## Calculate number of total oxygens, num epoxys and num hydroxyls 
## and return in dict
######################################################################
def manual_oxygen_num_type():
    OxygenNumAndType = {}
    # 30% hydroxyl; 70% dihydroxyls
    OxygenNumAndType["Hydroxyl"] = int(round(166*0.3,0))
    OxygenNumAndType["Dihydroxyl"] = int(round(166*0.7,0))
    if (OxygenNumAndType["Dihydroxyl"]%2) != 0:
        OxygenNumAndType["Dihydroxyl"] -= 1
        OxygenNumAndType["Hydroxyl"] += 1
    OxygenNumAndType["Epoxy"] = 129
    OxygenNumAndType["Total"] = 129+166

    print("Adding %s oxygen atoms:"%(OxygenNumAndType["Total"]))
    print("       O_hydroxyl:   %s"%(OxygenNumAndType["Hydroxyl"]))
    print("     O_dihydroxyl:   %s"%(OxygenNumAndType["Dihydroxyl"]))
    print("          O_epoxy:   %s"%(OxygenNumAndType["Epoxy"]))
    #print("    O_epoxy"%(NumO,NumHy,NumEp))

    return(OxygenNumAndType)

######################################################################
## Calculate number of total oxygens, num epoxys and num hydroxyls 
## and return in dict based on aging
## (Will only work if initial carbon lattice has a sufficient number 
##  of initial defect sites)
## TODO: Add warning if defect sites < edge sites
######################################################################
def calculate_oxygen_num_type_aged(NumC,ONumO,Temp,PercentO):
    # Carbon to oxygen ratio based on annealing paper (citation TBA)
    
    #if Temp == 600:
    #    PercentO = 0.87*PercentO
    #elif Temp == 900:
    #    PercentO = 0.71*PercentO
    #elif Temp == 1200:
    #    PercentO = 0.49*PercentO
    #elif Temp == 2500:
    #    PercentO = 0.14*PercentO

    if ONumO != 0:
        oxy = ONumO/NumC*100
        print("Reading structure with %s oxygen"%(oxy))
    NumO = int(NumC * (PercentO / 100 )) - ONumO
    OxygenNumAndType = {}

    # Built-in ratios based on annealing paper (citation TBA)
    ## Original percents: 21, 17, 4, 46, 5, 6, 0
    ## Modified to include carboxylic acid
    if Temp == 600:
        OxygenNumAndType["Dihydroxyl"] = int(round(NumO*0.19,0))
        OxygenNumAndType["Surface_Hydroxyl"] = int(round(NumO*0.16,0))
        OxygenNumAndType["Edge_Hydroxyl"] = int(round(NumO*0.04,0))
        OxygenNumAndType["Epoxy"] = int(round(NumO*0.43,0))
        OxygenNumAndType["Ether"] = int(round(NumO*0.05,0))
        OxygenNumAndType["Carbonyl"] = int(round(NumO*0.07,0))
        OxygenNumAndType["Carboxylic_Acid"] = int(round(NumO*0.06,0))  ## Percents not accurate at the moment
        #OxygenNumAndType["Carboxylic_Acid"] = 0  ## Percents not accurate at the moment
    ## Original percents: 17, 18, 5, 39, 10, 11, 0
    elif Temp == 900:
        OxygenNumAndType["Dihydroxyl"] = int(round(NumO*0.15,0))
        OxygenNumAndType["Surface_Hydroxyl"] = int(round(NumO*0.16,0))
        OxygenNumAndType["Edge_Hydroxyl"] = int(round(NumO*0.05,0))
        OxygenNumAndType["Epoxy"] = int(round(NumO*0.34,0))
        OxygenNumAndType["Ether"] = int(round(NumO*0.08,0))
        OxygenNumAndType["Carbonyl"] = int(round(NumO*0.11,0))
        OxygenNumAndType["Carboxylic_Acid"] = int(round(NumO*0.11,0))  ## Percents not accurate at the moment
        #OxygenNumAndType["Carboxylic_Acid"] = 0  ## Percents not accurate at the moment
    ## Original percents: 12, 14, 6, 22, 18, 29, 0
    elif Temp == 1200:
        OxygenNumAndType["Dihydroxyl"] = int(round(NumO*0.12,0))
        OxygenNumAndType["Surface_Hydroxyl"] = int(round(NumO*0.14,0))
        OxygenNumAndType["Edge_Hydroxyl"] = int(round(NumO*0.06,0))
        OxygenNumAndType["Epoxy"] = int(round(NumO*0.22,0))
        OxygenNumAndType["Ether"] = int(round(NumO*0.18,0))
        OxygenNumAndType["Carbonyl"] = int(round(NumO*0.14,0))
        OxygenNumAndType["Carboxylic_Acid"] = int(round(NumO*0.14,0))  ## Percents not accurate at the moment
        #OxygenNumAndType["Carboxylic_Acid"] = 0  ## Percents not accurate at the moment
    ## Original percents: 2, 2, 4, 15, 47, 31, 0
    elif Temp == 2500:
        OxygenNumAndType["Dihydroxyl"] = int(round(NumO*0.02,0))
        OxygenNumAndType["Surface_Hydroxyl"] = int(round(NumO*0.02,0))
        OxygenNumAndType["Edge_Hydroxyl"] = int(round(NumO*0.04,0))
        OxygenNumAndType["Epoxy"] = int(round(NumO*0.15,0))
        OxygenNumAndType["Ether"] = int(round(NumO*0.47,0))
        OxygenNumAndType["Carbonyl"] = int(round(NumO*0.15,0))
        OxygenNumAndType["Carboxylic_Acid"] = int(round(NumO*0.15,0))  ## Percents not accurate at the moment


    if (OxygenNumAndType["Dihydroxyl"]%2) != 0:
        OxygenNumAndType["Dihydroxyl"] -= 1
        OxygenNumAndType["Surface_Hydroxyl"] += 1


    OxygenNumAndType["Total"] = OxygenNumAndType["Dihydroxyl"] + OxygenNumAndType["Surface_Hydroxyl"] + OxygenNumAndType["Edge_Hydroxyl"] + OxygenNumAndType["Epoxy"] + OxygenNumAndType["Ether"] + OxygenNumAndType["Carbonyl"] + OxygenNumAndType["Carboxylic_Acid"]*2

    print("Adding %s oxygen atoms:"%(OxygenNumAndType["Total"]))
    print("       O_hydroxyl:   %s"%(OxygenNumAndType["Dihydroxyl"]))
    print("O_surfacehydroxyl:   %s"%(OxygenNumAndType["Surface_Hydroxyl"]))
    print("   O_edgehydroxyl:   %s"%(OxygenNumAndType["Edge_Hydroxyl"]))
    print("          O_epoxy:   %s"%(OxygenNumAndType["Epoxy"]))
    print("          O_ether:   %s"%(OxygenNumAndType["Ether"]))
    print("       O_carbonyl:   %s"%(OxygenNumAndType["Carbonyl"]))
    print("       O_carbacid:   %s"%(OxygenNumAndType["Carboxylic_Acid"]*2))

    return(OxygenNumAndType)

######################################################################
## Define defect edge sites
######################################################################
def define_defect_edge_sites(coords,atomtypes,abc):
    DefectEdgeSites = []
    NNList = nearest_neighbors_fv(coords, abc, 1.65)
    for i in range(0,len(coords)):
        if len(NNList[i]) < 3 and "C" in atomtypes[i]:
            DefectEdgeSites.append(i)
    return(DefectEdgeSites)

######################################################################
## All atoms are possible defect sites
######################################################################
def use_all_defect_sites(atomtypes):
    PossibleDefectSites = []
    for i in range(0,len(atomtypes)):
        if "C" in atomtypes[i]:
            PossibleDefectSites.append(i)
    return(PossibleDefectSites)

######################################################################
## All carbon atoms which do not have oxygen bound
######################################################################
def recalculate_all_defect_sites(coords,atomtypes,abc):
    PossibleDefectSites = []
    for i in range(0,len(atomtypes)):
        if "C" in atomtypes[i]:
            NNList = nearest_neighbor_single(coords[i],coords,abc,1.75)
            oxygen = False
            for atom in NNList:
                if "O" in atomtypes[atom]:
                    oxygen = True
            if oxygen != True:
                PossibleDefectSites.append(i)
    return(PossibleDefectSites)

######################################################################
## Add epoxy oxygen
## 2 C atoms where epoxy is added      o
## c1,c1  = [ x1, y1, z1 ],           / \
##          [ x2, y2, z2 ]          c1---c2
## Assumes c1 and c2 are in the xy plane
######################################################################
def add_epoxy(c1,c2,up,abc):
    if (c1[0] - c2[0]) >= (abc[0]*0.5):
        c1[0] -= abc[0]
    elif (c1[0] - c2[0]) <= -(abc[0]*0.5):
        c2[0] -= abc[0]
    if (c1[1] - c2[1]) >= (abc[0]*0.5):
        c1[1] -= abc[1]
    elif (c1[1] - c2[1]) <= -(abc[0]*0.5):
        c2[1] -= abc[1]
    ox = (c2[0]+c1[0])*0.5
    oy = (c2[1]+c1[1])*0.5
    oz = (c2[2]+c1[2])*0.5
    if up == True:
        oz += 1.283
    else:
        oz -= 1.283
    oCoords = [ox,oy,oz]
    return(oCoords)

######################################################################
## Add ether oxygen
## 2 C atoms where epoxy is added      o
## c1,c1  = [ x1, y1, z1 ],           / \
##          [ x2, y2, z2 ]          c1   c2    (o in xy plane)
## Assumes c1 and c2 are in the xy plane
######################################################################
def add_ether(c1,c2,pos,abc):
    oz = (c1[2] + c2[2])*0.5
    if (c1[0] - c2[0]) >= (abc[0]*0.5):
        c1[0] -= abc[0]
    elif (c1[0] - c2[0]) <= -(abc[0]*0.5):
        c2[0] -= abc[0]
    if (c1[1] - c2[1]) >= (abc[0]*0.5):
        c1[1] -= abc[1]
    elif (c1[1] - c2[1]) <= -(abc[0]*0.5):
        c2[1] -= abc[1]
    result = get_new_position(a=np.array([c1[0],c1[1]]),b=np.array([c2[0],c2[1]]),length=0.79)
    oCoords = [result['position_vect'][0], result['position_vect'][1], oz]
    midpoint = [result['midpoint_vect'][0], result['midpoint_vect'][1], oz]
    return(oCoords,midpoint)



######################################################################
## Apply rotation for bond
######################################################################
def rotate_atom(x,y,z):
    angle = np.pi/3
    possibleRotations = [0,angle,angle*2,angle*3,angle*4,angle*5]
    RandomRotation = random.choice(possibleRotations)
    x,y = rotate( (0,0), (x,y), RandomRotation)
    return(x,y,z)

######################################################################
## Add hydroxyl oxygen and hydrogen      c---o
## c = [ x1, y1, z1]                          \
##                                             h
## Assumes neighboring carbons are in xy plane
## OH bond is a bit funny to break any possible overlap
## C-O bond is LONG to prevent overlap - vmd likes to draw C-O bonds 
## apparently and for topology generation this is a problem
######################################################################
def add_hydroxyl(c,up):
    ox,oy,oz = c[0],c[1],c[2]
    if up == True:
        #oz += 1.32
        oz += 1.7
        #hx,hy,hz = rotate_atom(0.885,0.0,0.314)
        hx,hy,hz = rotate_atom(0.4,0.0,0.8)
    else:
        #oz -= 1.32
        oz -= 1.7
        #hx,hy,hz = rotate_atom(0.885,0.0,-0.314)
        hx,hy,hz = rotate_atom(0.4,0.0,-0.8)
    oCoords = [ox,oy,oz]
    hCoords = [ox+hx,oy+hy,oz+hz]
    return(oCoords,hCoords)

######################################################################
## Add dihydroxyl oxygen and hydrogen       h1     r
## c1 = [ x1, y1, z1]                         \    |
## c2 = [ x2, y2, z2]                         o1---c1
##                                                 |
##                                                 c2---o2
## Assumes neighboring carbons are in xy plane     |     \
##                                                 r      h2
######################################################################
def add_dihydroxyl(c1,c2,abc):
    if (c1[0] - c2[0]) >= (abc[0]*0.5):
        c1[0] -= abc[0]
    elif (c1[0] - c2[0]) <= -(abc[0]*0.5):
        c2[0] -= abc[0]
    if (c1[1] - c2[1]) >= (abc[0]*0.5):
        c1[1] -= abc[1]
    elif (c1[1] - c2[1]) <= -(abc[0]*0.5):
        c2[1] -= abc[1]
    o1x,o1y,o1z = c1[0],c1[1],c1[2]
    o2x,o2y,o2z = c2[0],c2[1],c2[2]
    up = coin_flip()
    if up == True:
        #o1z += 1.32
        #o2z -= 1.32
        o1z += 1.7
        o2z -= 1.7
        #h1x,h1y,h1z = rotate_atom(0.885,0.0,0.314)
        #h2x,h2y,h2z = rotate_atom(0.885,0.0,-0.314)
        h1x,h1y,h1z = rotate_atom(0.4,0.0,0.8)
        h2x,h2y,h2z = rotate_atom(0.4,0.0,-0.8)
    else:
        #o1z -= 1.32
        #o2z += 1.32
        o1z -= 1.7
        o2z += 1.7
        #h1x,h1y,h1z = rotate_atom(0.885,0.0,-0.314)
        #h2x,h2y,h2z = rotate_atom(0.885,0.0,0.314)
        h1x,h1y,h1z = rotate_atom(0.4,0.0,-0.8)
        h2x,h2y,h2z = rotate_atom(0.4,0.0,0.8)
        #h2x,h2y,h2z = rotate_atom(0.8054,0.4435,0.3091)
    o1Coords,o2Coords = [o1x,o1y,o1z],[o2x,o2y,o2z]
    h1Coords,h2Coords = [o1x+h1x,o1y+h1y,o1z+h1z],[o2x+h2x,o2y+h2y,o2z+h2z]
    return(o1Coords,h1Coords,o2Coords,h2Coords)

######################################################################
## Add carboxylic acid oxygen and hydrogen   h1---o1   o2
## c = [ x1, y1, z1]                                \ //
##                                               r---c1---r       
######################################################################
def add_carboxylic_acid(c):
    up = coin_flip()
    if up == True:
        o1x,o1y,o1z = rotate_atom(0.0,1.07,0.615)
        o2x,o2y,o2z = [-o1x,-o1y,o1z]
        h1x,h1y,h1z = rotate_atom(0.8054,0.4435,0.3091)
    else:
        o1x,o1y,o1z = rotate_atom(0.0,1.07,-0.615)
        o2x,o2y,o2z = [-o1x,-o1y,o1z]
        h1x,h1y,h1z = rotate_atom(0.8054,0.4435,-0.3091)
    o1Coords = [o1x+c[0],o1y+c[1],o1z+c[2]]
    o2Coords = [o2x+c[0],o2y+c[1],o2z+c[2]]
    h1Coords = [ox+hx,oy+hy,oz+hz]
    return(o1Coords,o2Coords,h1Coords)

######################################################################
## Add carbonyl oxygen        |    
## c = [ x1, y1, z1]          c=o
##                            |
######################################################################
def add_carbonyl(c,up):
    ox,oy,oz = c[0],c[1],c[2]
    if up == True:
        oz += 1.5
    else:
        oz -= 1.5
    oCoords = [ox,oy,oz]
    return(oCoords)

######################################################################
## Add carbonyl oxygen in the z-plane    \
## c = [ x1, y1, z1 ]                     c=o
##                                       /
## TODO: Will actually need other carbon info
######################################################################
#def add_carbonyl_inplane(c):
#    return(oCoords)

######################################################################
## Ensure that there are the same # of oxygens on either side
######################################################################
def choose_direction(o_up,o_down,o_tot):
    if o_up <= (o_tot*0.5) and o_down <= (o_tot*0.5):
        up = coin_flip()
    elif o_up >= (o_tot*0.5) and o_down <= (o_tot*0.5):
        up = False
    elif o_up >= (o_tot*0.5) and o_down <= (o_tot*0.5):
        up = True
    else: #Handle odd numbers
        up = coin_flip()
    if up == False:
        o_down += 1
    else:
        o_up += 1
    return(o_up,o_down,up)

def pick_epoxy(coords,NewO,NNList,OxygenNumAndType):
    for epoxy in range(0,OxygenNumAndType["Epoxy"]):
        c1index = random.choice(PossibleDefectSites)
        placed = False
        for neighbor in NNList[c1index]:
            if neighbor in PossibleDefectSites and placed == False:
                c2index = neighbor
                oCoords = add_epoxy(coords[c1index],coords[c2index])
                placed = True
                PossibleDefectSites.remove(c2index)
                PossibleDefectSites.remove(c2index)
        if placed == False:
            c2index = random.choice(NNList[c1index])
            oCoords = add_epoxy(coords[c1index],coords[c2index])
        PossibleDefectSites.remove(c1index)
        NewO.append(oCoords)
    return(PossibleDefectSites,NewO)

######################################################################
## Only include sites within a circle of a given radius from a designated
## site
######################################################################
def possible_sites_circle(coords,atomtypes,abc,site,radius):
    try:
        anchorx = coords[site][0]
        anchory = coords[site][1]
        anchorz = coords[site][2]
    ## use specific coordinates
    except IndexError:
        print(site)
        anchorx = float(site[0])
        anchory = float(site[1])
        anchorz = float(site[2])
    anchor = [anchorx,anchory,anchorz]
    DefectPossibleNums = []

    for i in range(0,len(coords)):
        dx = abs(anchorx - coords[i][0])
        if dx >= (abc[0]*0.5):
            dx -= abc[0]
        dy = abs(anchory - coords[i][1])
        if dy >= (abc[1]*0.5):
            dy -= abc[1]
        dist = (dx**2+dy**2)**0.5
        if dist <= radius:
            DefectPossibleNums.append(i)
    return(DefectPossibleNums)

######################################################################
## Check that carbons have 1 NN within PossibleDefectSites
######################################################################
def assess_two_carbon_group(coords,PossibleDefectSites,NNList):
    twocarbongroupsites = []
    for atom in PossibleDefectSites:
        for neighbor in NNList[atom]:
            if neighbor in PossibleDefectSites:
                if atom not in twocarbongroupsites:
                    twocarbongroupsites.append(atom)
    return(twocarbongroupsites)

######################################################################
## Remove triple coordinated O by checking if neighboring carbons 
## are also neighboring carbons
######################################################################
def remove_triple_coordinated_o(EtherSites,EtherNNList):
    #for Site in EtherSites:
    #    print(Site,end=" ")
    #print("")
    NewEtherSites = []
    for Site in EtherSites:
        NewEtherSites.append(Site)
    for Site in EtherSites:
        NumNeighbors = 0
        for Neighbor in EtherNNList[Site]:
            for SecondNeighbor in EtherNNList[Neighbor]:
                if SecondNeighbor in EtherNNList[Site]:
                    NumNeighbors += 1
        if NumNeighbors >= 2:
            NewEtherSites.remove(Site)
    #for Site in NewEtherSites:
    #    print(Site,end=" ")
    #print("")
    return(NewEtherSites)

######################################################################
## Ether sites are defect atoms with < 3 NN 
## and which have at least 1 atom (that also has <3 NN) within 2.9 Ang
## To remove triple coordinated oxygen, defect atoms which have the 
## same third atom are removed
## exclude pbc boundaries as these cause issues
######################################################################
def find_ether_sites(coords,abc,EdgeSites,NNList):
    EtherPossibleSites = []
    for atom in EdgeSites:
        ix = coords[atom][0]
        iy = coords[atom][1]
        iz = coords[atom][2]
        if (ix < 1.0) or (ix > (abc[0]-1.0)) or (iy < 1.0) or (iy > (abc[1]-1.0)):
            pass
        elif len(NNList[atom]) < 3:
            EtherPossibleSites.append(atom)
    EtherNNList = {}
    for i in EtherPossibleSites:
        ix = coords[i][0]
        iy = coords[i][1]
        iz = coords[i][2]
        NeighborElements = []
        for j in EtherPossibleSites:
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
            if d <= 2.7 and d >= 1.7 and i != j:
                NeighborElements.append(j)
        EtherNNList[i] = NeighborElements #Dict because there are a different number of atoms 
    EtherSites = []
    CarbAcidInverseSites = []
    for i in EtherPossibleSites:
        if len(EtherNNList[i]) != 0:
            EtherSites.append(i)
            CarbAcidInverseSites.append(i)
    EtherSites = remove_triple_coordinated_o(EtherSites,EtherNNList)
    #for i in EtherSites:
        #print("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",250,"N","UNL",1,coords[i][0],coords[i][1],coords[i][2]+31.675,1.00,0.00,"N"))
    return(EtherSites,EtherNNList,CarbAcidInverseSites)

######################################################################
## Check if ether O overlaps with carbon
######################################################################
def check_ether_overlap(oCoords,c1index,c2index,coords,abc,NNList,NewO):
    #Check if oxygen overlaps with NN on 1st carbon
    for neighbor in NNList[c1index]:
        dx = abs(oCoords[0] - coords[neighbor][0])
        if dx >= (abc[0]*0.5):
            dx -= abc[0]
        dy = abs(oCoords[1] - coords[neighbor][1])
        if dy >= (abc[0]*0.5):
            dy -= abc[0]
        dz = abs(oCoords[2] - coords[neighbor][2])
        if dz >= (abc[0]*0.5):
            dz -= abc[0]
        d = ( dx**2 + dy**2 + dz**2)**0.5
        if d < 1.9:
            return(True)
    #Check if oxygen overlaps with NN on 2nd carbon
    for neighbor in NNList[c2index]:
        dx = abs(oCoords[0] - coords[neighbor][0])
        if dx >= (abc[0]*0.5):
            dx -= abc[0]
        dy = abs(oCoords[1] - coords[neighbor][1])
        if dy >= (abc[0]*0.5):
            dy -= abc[0]
        dz = abs(oCoords[2] - coords[neighbor][2])
        if dz >= (abc[0]*0.5):
            dz -= abc[0]
        d = ( dx**2 + dy**2 + dz**2)**0.5
        if d < 1.9:
            return(True)
    for oxygen in NewO:
        dx = abs(oCoords[0] - oxygen[0])
        if dx >= (abc[0]*0.5):
            dx -= abc[0]
        dy = abs(oCoords[1] - oxygen[1])
        if dy >= (abc[0]*0.5):
            dy -= abc[0]
        dz = abs(oCoords[2] - oxygen[2])
        if dz >= (abc[0]*0.5):
            dz -= abc[0]
        d = ( dx**2 + dy**2 + dz**2)**0.5
        if d < 1.9:
            return(True)
    return(False)

######################################################################
## Add all ether
######################################################################
def add_all_ether(NewO,OxygenNumAndType,coords,abc,EdgeSites,NNList):
    EtherSites,EtherNNList,CarbAcidInverseSites = find_ether_sites(coords,abc,EdgeSites,NNList)
    CarbAcidSites = []
    for i in EdgeSites:
        if i not in CarbAcidInverseSites:
            CarbAcidSites.append(i)
    NonPlacedEther = 0
    for ether in range(0,OxygenNumAndType["Ether"]):
        placed = False
        count = 0
        while placed == False and count < 100:
            count += 1
            try:
                c1index = random.choice(EtherSites)
                #print(c1index)
            except IndexError:
                print("Warning: there are no more carbons available for 2 site ether placements. Ending ether placement early (%s placed out of %s)."%((ether-NonPlacedEther),OxygenNumAndType["Ether"]))
                return(NewO,EdgeSites,CarbAcidSites)
            placed = False
            for neighbor in EtherNNList[c1index]:
                if neighbor in EtherSites and placed == False:
                    c2index = neighbor
                    pos = random.choice([True,False])
                    oCoords,midpoint = add_ether(coords[c1index],coords[c2index],pos,abc)
                    EtherOverlap = check_ether_overlap(oCoords,c1index,c2index,coords,abc,NNList,NewO)
                    if EtherOverlap == False:
                        placed = True
                        NewO.append(oCoords)
                    elif EtherOverlap == True:
                        if pos == True:
                            pos = False
                        else:
                            pos = True
                        oCoords_xy = rotate( (midpoint[0],midpoint[1]), (oCoords[0],oCoords[1]), np.pi)
                        oCoords = [ oCoords_xy[0], oCoords_xy[1], oCoords[2] ]
                        EtherOverlap = check_ether_overlap(oCoords,c1index,c2index,coords,abc,NNList,NewO)
                        if EtherOverlap == False:
                            placed = True
                            NewO.append(oCoords)
                        elif EtherOverlap == True:
                            try:
                                EtherSites.remove(c1index)  ## Removes problem atoms at the cost of over-removing
                                CarbAcidSites.append(c1index)
                            except ValueError:
                                pass
                            try:
                                EtherSites.remove(c2index)
                                CarbAcidSites.append(c2index)
                            except ValueError:
                                pass
                            #EtherNNList[c1index] = EtherNNList[c1index].remove(neighbor)
        if count == 100:
            NonPlacedEther += 1
        try:
            EtherSites.remove(c1index)
        except ValueError:
            pass
        try:
            EtherSites.remove(c2index)
        except ValueError:
            pass
        try:
            EdgeSites.remove(c1index)
        except ValueError:
            pass
        try:
            EdgeSites.remove(c2index)
        except ValueError:
            pass
        #NewO.append(oCoords)
        ## remove overlapping sites from carb acid placement
        for neighbor in nearest_neighbor_single(oCoords,coords,abc,2.5):
            try:
                CarbAcidSites.remove(neighbor)
            except ValueError:
                 pass
    if NonPlacedEther > 0:
            print("Warning: not all ethers placed (%s not placed out of %s)"%(NonPlacedEther,OxygenNumAndType["Ether"]))
    for i in EtherSites:
        CarbAcidSites.append(i)
    return(NewO,EdgeSites,CarbAcidSites)

######################################################################
## Add all carbonyl
######################################################################
def add_all_carbonyl(NewO,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down):
    for carbonyl in range(0,OxygenNumAndType["Carbonyl"]):
        try:
            c1index = random.choice(EdgeSites)
        except IndexError:
            print("Warning: All possible edge defect sites have been used.")
            return(NewO,EdgeSites)
        o_up,o_down,up = choose_direction(o_up,o_down,OxygenNumAndType["Total"])
        #o1Coords = add_carbonyl_inplane(coords[c1index],up)
        o1Coords = add_carbonyl(coords[c1index],up)
        EdgeSites.remove(c1index)
        #Remove neighboring atoms from selection list
        for neighbor in NNList[c1index]:
            try:
                EdgeSites.remove(neighbor)
            except ValueError:
                pass
        #print(c1index)
        NewO.append(o1Coords)
    return(NewO,EdgeSites,o_up,o_down)

######################################################################
## Add all edge hydroxyls
######################################################################
def add_all_edgehydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down):
    for hydroxyl in range(0,OxygenNumAndType["Edge_Hydroxyl"]):
        placed = False
        while placed == False:
            try:
                c1index = random.choice(EdgeSites)
            except IndexError:
                print("Warning: All possible edge defect sites have been used. Placed %s out of %s edge hydroxyls."%(hydroxyl,OxygenNumAndType["Edge_Hydroxyl"]))
                return(NewO,NewH,EdgeSites,o_up,o_down)
            #o_up,o_down,up = choose_direction(o_up,o_down,OxygenNumAndType["Total"])
            #o1Coords,h1Coords = add_hydroxyl(coords[c1index],up)
            o_up,o_down,c1index,o1Coords,h1Coords,placed = attempt_hydroxyl_placement(c1index,NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
            EdgeSites.remove(c1index)
            #print(c1index)
        NewO.append(o1Coords)
        NewH.append(h1Coords)
    return(NewO,NewH,EdgeSites,o_up,o_down)

######################################################################
## Basically the same as the hydrogen placement function
######################################################################
def add_carboxylicacid_carbon(c1,c2,c3):
    # move to origin
    ba = np.array(c2) - np.array(c1)
    bc = np.array(c3) - np.array(c1)
    new = np.mean( np.array([ ba, bc ]), axis=0 )
    new = new*np.array([-1,-1])
    theta = abs(np.arctan((new[1]/new[0])))
    new_scaled = np.array([ np.cos(theta)*1.3, np.sin(theta)*1.3 ]) #for bond length of 1.3
    if new[0] < 0:
        new_scaled[0] = new_scaled[0]*-1
    if new[1] < 0:
        new_scaled[1] = new_scaled[1]*-1
    new_scaled += np.array(c1)
    return(new_scaled[0],new_scaled[1])

######################################################################
## Calculate coordinates for carbon of carb acid group
######################################################################
def choose_carbacid_atoms(c1index,coords,ABC,NNList):
    c1x,c1y,c1z = coords[c1index]
    c2 = [ coords[NNList[c1index][0]][0], coords[NNList[c1index][0]][1] ]
    c3 = [ coords[NNList[c1index][1]][0], coords[NNList[c1index][1]][1] ]
    c2,c3 = remove_PBC([c1x,c1y],c2,c3,ABC)
    c4x,c4y = add_carboxylicacid_carbon([c1x,c1y],c2,c3 )
    c4Coords = [c4x,c4y,c1z]
    return(c4Coords)

######################################################################
## Check that new carb acid carbon is only bonded to one other carbon
######################################################################
def check_carbacid_overlap(c1index,c4coords,ABC,coords,NewC,NewO,CarbAcidSites):
    NNList = nearest_neighbor_single(c4coords,coords,ABC,1.9)
    NNListNewC = nearest_neighbor_single(c4coords,NewC,ABC,1.9)
    NNListNewO = nearest_neighbor_single(c4coords,NewO,ABC,1.9)
    if (len(NNList)+len(NNListNewC)+len(NNListNewO)) > 1:
        CarbAcidSites.remove(c1index)
        return(True,CarbAcidSites)
    return(False,CarbAcidSites)

######################################################################
## Add all carboxylic acids
######################################################################
def add_all_carboxylicacid(ABC,NewC,NewO,NewH,OxygenNumAndType,coords,EdgeSites,CarbAcidSites,NNList):
    for carbox in range(0,OxygenNumAndType["Carboxylic_Acid"]):
        try:
            c1index = random.choice(CarbAcidSites)
        except IndexError:
            print("Warning: All possible carb acid sites have been used.")
            return(NewC,NewO,NewH,EdgeSites)

        c4Coords = choose_carbacid_atoms(c1index,coords,ABC,NNList)
        overlap,CarbAcidSites = check_carbacid_overlap(c1index,c4Coords,ABC,coords,NewC,NewO,CarbAcidSites) 

        while overlap == True:
            try:
                c1index = random.choice(CarbAcidSites)
            except IndexError:
                print("Warning: All possible carb acid sites have been used.")
                return(NewC,NewO,NewH,EdgeSites)

            c4Coords = choose_carbacid_atoms(c1index,coords,ABC,NNList)
            overlap,CarbAcidSites = check_carbacid_overlap(c1index,c4Coords,ABC,coords,NewC,NewO,CarbAcidSites)

        NewC.append(c4Coords)
        if coin_flip() == True:
            o1Coords = add_carbonyl(c4Coords,True)
            o2Coords,h2Coords = add_hydroxyl(c4Coords,False)
        else:
            o1Coords = add_carbonyl(c4Coords,False)
            o2Coords,h2Coords = add_hydroxyl(c4Coords,True)
        EdgeSites.remove(c1index)
        CarbAcidSites.remove(c1index)
        #Remove neighboring atoms from selection list
        for neighbor in NNList[c1index]:
            try:
                EdgeSites.remove(neighbor)
                CarbAcidSites.remove(neighbor)
            except ValueError:
                pass
        #print(c1index)
        for neighbor in nearest_neighbor_single(coords[c1index],coords,ABC,2.9):
            try:
                CarbAcidSites.remove(neighbor)
            except ValueError:
                pass
        NewO.append(o1Coords)
        NewO.append(o2Coords)
        NewH.append(h2Coords)
    return(NewC,NewO,NewH,EdgeSites)

######################################################################
## Add all dihydroxyls
######################################################################
def add_all_dihydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down):
    half = int(OxygenNumAndType["Dihydroxyl"]*0.5)
    twocarbongroupsites = assess_two_carbon_group(coords,PossibleDefectSites,NNList)
    for dihydroxyl in range(0,half):
        #twocarbongroupsites = assess_two_carbon_group(coords,PossibleDefectSites,NNList)
        c1index = random.choice(twocarbongroupsites)
        placed = False
        for neighbor in NNList[c1index]:
            if neighbor in twocarbongroupsites and placed == False:
                c2index = neighbor
                o1Coords,h1Coords,o2Coords,h2Coords = add_dihydroxyl(coords[c1index],coords[c2index],abc)
                placed = True
                #PossibleDefectSites.remove(c2index)
        #c1index = random.choice(PossibleDefectSites)
        if placed == False:
            c2index = random.choice(NNList[c1index])
            o1Coords,h1Coords,o2Coords,h2Coords = add_dihydroxyl(coords[c1index],coords[c2index],abc)
            print("Warning: overlapping atoms may occur (di-hydroxyl placement failure)")
        PossibleDefectSites.remove(c1index)
        PossibleDefectSites.remove(c2index)
        #Remove neighboring atoms from selection list
        for neighbor in NNList[c1index]:
            try:
                PossibleDefectSites.remove(neighbor)
            except ValueError:
                pass
        for neighbor in NNList[c2index]:
            try:
                PossibleDefectSites.remove(neighbor)
            except ValueError:
                pass
        for neighbor in NNList[c1index]:
            try:
                twocarbongroupsites.remove(neighbor)
            except ValueError:
                pass
        for neighbor in NNList[c2index]:
            try:
                twocarbongroupsites.remove(neighbor)
            except ValueError:
                pass
        NewO.append(o1Coords)
        NewO.append(o2Coords)
        NewH.append(h1Coords)
        NewH.append(h2Coords)
        o_up += 1
        o_down += 1
    return(NewO,NewH,PossibleDefectSites,o_up,o_down)

######################################################################
## Add all epoxys
######################################################################
def add_all_epoxy(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,NNListSW,o_up,o_down):
    twocarbongroupsites = assess_two_carbon_group(coords,PossibleDefectSites,NNList)
    for epoxy in range(0,OxygenNumAndType["Epoxy"]):
        #twocarbongroupsites = assess_two_carbon_group(coords,PossibleDefectSites,NNList)
        try:
            c1index = random.choice(twocarbongroupsites)
        except IndexError:
            print("There are no more carbons available for 2 site placements. Placing remaining epoxy groups (%s/%s) as hydroxyl."%(OxygenNumAndType["Epoxy"]-epoxy,OxygenNumAndType["Epoxy"]))
            return(NewO,NewH,PossibleDefectSites,o_up,o_down,OxygenNumAndType["Epoxy"]-epoxy)
            #c1index = random.choice(PossibleDefectSites)
        o_up,o_down,up = choose_direction(o_up,o_down,OxygenNumAndType["Total"])
        placed = False
        for neighbor in NNList[c1index]:
            if neighbor in twocarbongroupsites and placed == False:
                c2index = neighbor
                oCoords = add_epoxy(coords[c1index],coords[c2index],up,abc)
                placed = True
        if placed == False:
            c2index = random.choice(NNList[c1index])
            oCoords = add_epoxy(coords[c1index],coords[c2index],up,abc)
            print("Warning: overlapping atoms may occur (epoxy placement failure for epoxy #%s)"%(epoxy))
        #print(c1index,end=" ")
        #print(c2index,end=" ")
        #PossibleDefectSites.remove(c1index)
        #PossibleDefectSites.remove(c2index)
        #Remove neighboring atoms from selection list - allow neighboring atoms for epoxy
        for neighbor in NNList[c1index]:
            try:
                PossibleDefectSites.remove(neighbor)
            except ValueError:
                pass
        for neighbor in NNList[c2index]:
            try:
                PossibleDefectSites.remove(neighbor)
            except ValueError:
                pass
        for neighbor in NNList[c1index]:
            try:
                twocarbongroupsites.remove(neighbor)
            except ValueError:
                pass
        for neighbor in NNList[c2index]:
            try:
                twocarbongroupsites.remove(neighbor)
            except ValueError:
                pass
        NewO.append(oCoords)
    #print("")
    return(NewO,NewH,PossibleDefectSites,o_up,o_down,0)


######################################################################
## Add all hydroxyls
######################################################################
def add_all_hydroxyl(NewO,NewH,OxygenNumAndType,coords,PossibleDefectSites,NNList,o_up,o_down):
    warning_count = 0
    #UsedCarbons = []
    for hydroxyl in range(0,OxygenNumAndType["Surface_Hydroxyl"]):
        warning = False
        try:
            c1index = random.choice(PossibleDefectSites)
        except IndexError:
            #print("Warning: All possible defect sites have been used. Overlaps may occur.")
            warning_count += 1
            warning = True
        if warning == False:
            o_up,o_down,up = choose_direction(o_up,o_down,OxygenNumAndType["Total"])
            o1Coords,h1Coords = add_hydroxyl(coords[c1index],up)
            #UsedCarbons.append(c1index)
            PossibleDefectSites.remove(c1index)
            #Remove neighboring atoms from selection list
            for neighbor in NNList[c1index]:
                try:
                    PossibleDefectSites.remove(neighbor)
                except ValueError:
                    pass
            #print(c1index)
            NewO.append(o1Coords)
            NewH.append(h1Coords)
        #else:
            #PossibleDefectSites = recalculate_defect_sites(NewO,NewH,coords)
    #if warning_count != 0:
        #print("Warning: All possible defect sites have been used. Failure for %s placements."%(warning_count))
        #return(NewO,NewH,UsedCarbons,o_up,o_down)
    return(NewO,NewH,PossibleDefectSites,o_up,o_down)

######################################################################
## 
######################################################################
def check_overlap(o1Coords,h1Coords,NewO,NewH,coords,abc):
    ## just check oxygen for now
    ix,iy,iz = o1Coords
    NeighborElements = []
    for j in range(0,len(NewO)):
        jx = NewO[j][0]
        jy = NewO[j][1]
        jz = NewO[j][2]
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
        if d <= 1.75:
            NeighborElements.append(j)
            return(True)
    return(False)

######################################################################
## Place the hydroxyl, check if it overlaps with already placed O/H 
## atoms, if it does try the other side
######################################################################
def attempt_hydroxyl_placement(c1index,NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down):
    o_up_new,o_down_new,up = choose_direction(o_up,o_down,OxygenNumAndType["Total"])
    o1Coords,h1Coords = add_hydroxyl(coords[c1index],up)
    ## first try, placed with no changes
    if (check_overlap(o1Coords,h1Coords,NewO,NewH,coords,abc)) == False:
        return(o_up_new,o_down_new,c1index,o1Coords,h1Coords,True)
    ## try switching direction ##TODO: check that counters are counting direction properly
    if up == True:
        up = False
        o_down_new = o_down + 1
        o_up_new = o_up
    else:
        up = True
        o_up_new = o_up + 1
        o_down_new = o_down
    o1Coords,h1Coords = add_hydroxyl(coords[c1index],up)
    if (check_overlap(o1Coords,h1Coords,NewO,NewH,coords,abc)) == False:
        return(o_up_new,o_down_new,c1index,o1Coords,h1Coords,True)
    return(o_up,o_down,c1index,o1Coords,h1Coords,False)

def recalculate_defect_sites(NewO,coords,abc):
    PossibleDefectSites = []
    for i in range(0,len(coords)):
        ix = coords[i][0]
        iy = coords[i][1]
        iz = coords[i][2]
        NeighborElements = []
        for j in range(0,len(NewO)):
            jx = NewO[j][0]
            jy = NewO[j][1]
            jz = NewO[j][2]
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
            if d <= 1.7 and i != j:
                NeighborElements.append(j)
        #NNList.append(NeighborElements)
        if len(NeighborElements) == 0:
            PossibleDefectSites.append(i)
    return(PossibleDefectSites)

def recalculate_defect_sites_excluding_edges(PossibleDefectSites,coords,atomtypes,abc):
    NewPossibleDefectSites = []
    for i in PossibleDefectSites:
        if (coords[i][0] > 1.0) and (coords[i][0] < (abc[0]-1.0)) and (coords[i][1] > 1.0) and (coords[i][1] < (abc[1] - 1.0)):
            NewPossibleDefectSites.append(i)
    return(NewPossibleDefectSites)

######################################################################
## Add all hydroxyls
## Check nearest neighbors for every placement 
## (slower, but allows more oxygen placment)
######################################################################
def add_all_hydroxyl_unrestricted(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down):
    ## Recalculate possible defect sites as all carbons that don't have an oxygen bound
    PossibleDefectSites = recalculate_defect_sites(NewO,coords,abc)

    for hydroxyl in range(0,OxygenNumAndType["Surface_Hydroxyl"]):
        placed = False
        while placed == False:
            try:
                c1index = random.choice(PossibleDefectSites)
            except IndexError:
                print("Warning: All possible defect sites have been used. Exiting placement.")
                return(NewO,NewH,PossibleDefectSites,o_up,o_down)
            o_up,o_down,c1index,o1Coords,h1Coords,placed = attempt_hydroxyl_placement(c1index,NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down)
            PossibleDefectSites.remove(c1index)
        NewO.append(o1Coords)
        NewH.append(h1Coords)
    return(NewO,NewH,PossibleDefectSites,o_up,o_down)

######################################################################
## Add new atoms
######################################################################
def add_new_atoms(coords,oldatomtypes,NewC,NewO,NewH):
    newcoords = coords
    atomtypes = []
    for oldatom in range(0,len(coords)):
    #    print("C %s %s %s"%(cCoords[0],cCoords[1],cCoords[2]))
         #newcoords.append(coords[oldatom])
         atomtypes.append(oldatomtypes[oldatom])
    for cCoords in NewC:
    #    print("O %s %s %s"%(oCoords[0],oCoords[1],oCoords[2]))
         newcoords = np.vstack([newcoords, cCoords])
         atomtypes.append("C")
    for oCoords in NewO:
    #    print("O %s %s %s"%(oCoords[0],oCoords[1],oCoords[2]))
         newcoords = np.vstack([newcoords, oCoords])
         atomtypes.append("O")
    for hCoords in NewH:
    #    print("H %s %s %s"%(hCoords[0],hCoords[1],hCoords[2]))
         newcoords = np.vstack([newcoords, hCoords])
         atomtypes.append("H")
    return(newcoords,atomtypes)

def remove_PBC(c1,c2,c3,abc):
    #print(c1,c2,c3)
    ## Center c2 and c3 about c1
    d3 = [ c3[0] - c1[0], c3[1] - c1[1] ]
    if d3[0] >= (abc[0]*0.5):
        c3[0] -= abc[0]
    elif d3[0] <= -(abc[0]*0.5):
        c3[0] += abc[0]
    if d3[1] >= (abc[1]*0.5):
        c3[1] -= abc[1]
    elif d3[1] <= -(abc[1]*0.5):
        c3[1] += abc[1]
    d2 = [ c2[0] - c1[0], c2[1] - c1[1] ]
    if d2[0] >= (abc[0]*0.5):
        c2[0] -= abc[0]
    elif d2[0] <= -(abc[0]*0.5):
        c2[0] += abc[0]
    if d2[1] >= (abc[1]*0.5):
        c2[1] -= abc[1]
    elif d2[1] <= -(abc[1]*0.5):
        c2[1] += abc[1]
    #print(c1,c2,c3)
    return(c2,c3)

######################################################################
## Average the two carbon positions, flip it across the original carbon 
## and scale the bond to 1.0
######################################################################
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
    if new[0] == 0.0:
        theta = np.pi*0.5
    else:
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
    NNList = nearest_neighbors_fv(Coords,ABC,1.75)
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
            c2,c3 = remove_PBC([c1x,c1y],c2,c3,ABC)
            ## Calculate hydrogen coordinates
            hx,hy = calculate_hydrogen_coordinates([c1x,c1y],c2,c3)
            hCoords = [ hx,hy,Coords[Atom][2] ]
            ## Append coordinates
            NewH.append(hCoords)
    return(NewH)

######################################################################
## Wrap coordinates back to within the periodic boundary conditions
## (ignore z)
######################################################################
def wrap_coords(coords,abc):
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
def recalculate_edge_sites(Coords,ABC,atomtypes,includeNN):
    NNList = nearest_neighbors_fv(Coords,ABC,1.75)
    EdgeSites = []
    print(len(Coords))
    print(len(atomtypes))
    for Atom in range(0,len(Coords)):
        Neighbors = NNList[Atom]
        CNeighbors = []
        #print(Neighbors)
        for Neighbor in Neighbors:
            #print(Neighbor)
            if "C" in atomtypes[Neighbor]:
                CNeighbors.append(Neighbor)
        if len(CNeighbors) < 2 and "C" in atomtypes[Atom]:
            #print("< 2 for %s"%(Atom))
            pass
        elif len(Neighbors) < 3 and "C" in atomtypes[Atom]:
            EdgeSites.append(Atom)

    if includeNN == True:
        #for i in EdgeSites:
        #    print(i,end=" ")
        NewEdgeSites = []
        for i in EdgeSites:
            NewEdgeSites.append(i)
        for i in EdgeSites:
            for nn in NNList[i]:
                if (nn not in EdgeSites) and (nn not in NewEdgeSites):
                    bound_oxygen = False
                    carb_acid = False
                    for nnn in NNList[nn]:
                        if nnn != i:
                            if "O" in atomtypes[nnn]:
                                bound_oxygen = True
                            ## exclude carb acid carbon
                            for nnnn in NNList[nnn]:
                                if nnnn != nn:
                                    if "C" not in atomtypes[nnnn]:
                                        carb_acid = True
                    if bound_oxygen == False and carb_acid == False:
                        NewEdgeSites.append(nn)
        return(NewEdgeSites)
    return(EdgeSites)

######################################################################
## Add oxygens and hydrogens
######################################################################
def add_oxygen_hydrogen(coords,atomtypes,abc,PossibleDefectSites,OxygenNumAndType):
    NNList = nearest_neighbors_fv(coords,abc,1.75)
    ## SW
    NNListSW = nearest_neighbors_fv(coords,abc,1.4)
    #if aged == False:
    NewO = []
    NewC = []
    NewH = []
    # Counters for O direction
    o_up = 0
    o_down = 0
    if len(PossibleDefectSites) >= OxygenNumAndType["Total"]:
        NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_dihydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down)
        NewO,NewH,PossibleDefectSites,o_up,o_down,unplacedepoxy = add_all_epoxy(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,NNListSW,o_up,o_down)
        OxygenNumAndType["Surface_Hydroxyl"] += unplacedepoxy
        NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_hydroxyl(NewO,NewH,OxygenNumAndType,coords,PossibleDefectSites,NNList,o_up,o_down)
        OxygenNumAndType["Surface_Hydroxyl"] = OxygenNumAndType["Total"] - len(NewO)
        NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_hydroxyl_unrestricted(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down)
    else:
        print("less possible defect sites than total oxygen to be added. functionality to be added")
    #print(len(coords)+len(NewO)+len(NewH))
    newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
    NewH = add_hydrogen(newcoords,newatomtypes,abc,NewH)
    newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
    newcoords = wrap_coords(newcoords,abc)
    return(newcoords,newatomtypes)

######################################################################
## Add oxygens and hydrogens (aged)
######################################################################
def add_oxygen_hydrogen_aged(coords,atomtypes,abc,PossibleDefectSites,EdgeSites,OxygenNumAndType):
    ## Remove edge sites from possibledefectsites
    for Site in EdgeSites:
        if Site in PossibleDefectSites:
            PossibleDefectSites.remove(Site)
    NNList = nearest_neighbors_fv(coords,abc,1.75)
    NNListSW = nearest_neighbors_fv(coords,abc,1.4)
    NewC = []
    NewO = []
    NewH = []
    # Counters for O direction
    o_up = 0
    o_down = 0
    if len(PossibleDefectSites) >= OxygenNumAndType["Total"]:
        NewO,EdgeSites,CarbAcidSites = add_all_ether(NewO,OxygenNumAndType,coords,abc,EdgeSites,NNList)
        #for i in CarbAcidSites:
            #print(i,end=" ")
        #write_pdb_file(coords,atomtypes,abc,"ether.pdb")
        #newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
        #write_pdb_file(newcoords,atomtypes,abc,"ether.pdb")
        NewC,NewO,NewH,EdgeSites = add_all_carboxylicacid(abc,NewC,NewO,NewH,OxygenNumAndType,coords,EdgeSites,CarbAcidSites,NNList)
        NewO,EdgeSites,o_up,o_down = add_all_carbonyl(NewO,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
        NewO,NewH,EdgeSites,o_up,o_down = add_all_edgehydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
        NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_dihydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down)
        NewO,NewH,PossibleDefectSites,o_up,o_down,unplacedepoxy = add_all_epoxy(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,NNListSW,o_up,o_down)
        OxygenNumAndType["Surface_Hydroxyl"] += unplacedepoxy
        NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_hydroxyl(NewO,NewH,OxygenNumAndType,coords,PossibleDefectSites,NNList,o_up,o_down)
    newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
    NewH = add_hydrogen(newcoords,newatomtypes,abc,NewH)
    newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
    newcoords = wrap_coords(newcoords,abc)
    return(newcoords,newatomtypes)

######################################################################
## Add oxygens and hydrogens (aged) where all edge sites are used 
## first
######################################################################
def add_oxygen_hydrogen_aged_edge(coords,atomtypes,abc,PossibleDefectSites,EdgeSites,OxygenNumAndType):
    ## Remove edge sites from possibledefectsites
    for Site in EdgeSites:
        if Site in PossibleDefectSites:
            PossibleDefectSites.remove(Site)
    NNList = nearest_neighbors_fv(coords,abc,1.75)
    NNListSW = nearest_neighbors_fv(coords,abc,1.4)
    NewC = []
    NewO = []
    NewH = []
    # Counters for O direction
    o_up = 0
    o_down = 0
    if len(PossibleDefectSites) >= OxygenNumAndType["Total"]:
        ## Add ether to edge sites
        NewO,EdgeSites,CarbAcidSites = add_all_ether(NewO,OxygenNumAndType,coords,abc,EdgeSites,NNList)
        NumEther = len(NewO)
        print("Added %s Ether oxygen groups"%(NumEther))

        ## Add carb acid to edge sites
        NewC,NewO,NewH,EdgeSites = add_all_carboxylicacid(abc,NewC,NewO,NewH,OxygenNumAndType,coords,EdgeSites,CarbAcidSites,NNList)
        NumCarbAcid = len(NewO) - NumEther
        print("Added %s Carb Acid groups"%(int(NumCarbAcid*0.5)))

        ## Add carbonyl to edge sites
        NewO,EdgeSites,o_up,o_down = add_all_carbonyl(NewO,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
        NumCarbonyl = len(NewO) - NumEther - NumCarbAcid
        print("Added %s Carbonyl groups"%(NumCarbonyl))

        ## Recalculate edge sites to include NN
        newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
        EdgeSites = recalculate_edge_sites(newcoords,abc,newatomtypes,True)
        write_pdb_file(newcoords,newatomtypes,abc,"edgecheck.pdb")
        print("Updated Edge Sites:")
        for Site in EdgeSites:
            print(Site,end=" ")
            try:
                PossibleDefectSites.remove(Site)
            except ValueError:
                pass
        print("")

        ## Add epoxy to edge sites
        NewO,NewH,EdgeSites,o_up,o_down,unplacedepoxy = add_all_epoxy(NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,NNListSW,o_up,o_down)
        NumEdgeEpoxy = len(NewO) - NumEther - NumCarbAcid - NumCarbonyl
        print("Added %s Edge Epoxy groups"%(NumEdgeEpoxy))
        if unplacedepoxy > 0:
            ## Add remaining epoxy to surface sites
            NewO,NewH,PossibleDefectSites,o_up,o_down,unplacedepoxy = add_all_epoxy(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,NNListSW,o_up,o_down)
            NumSurfaceEpoxy = len(NewO) - NumEther - NumCarbAcid - NumCarbonyl - NumEdgeEpoxy
            print("Added %s Surface Epoxy groups"%(NumSurfaceEpoxy))
            #OxygenNumAndType["Edge_Hydroxyl"] += OxygenNumAndType["Dihydroxyl"]
            #NumEdgeDiHydrox = 0
        else:
            ## Add dihydroxyl to edge sites
            #NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_dihydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
            #NumEdgeDiHydrox = len(NewO) - NumEther - NumCarbAcid - NumCarbonyl - NumEdgeEpoxy
            #print("Added %s Edge DiHydrox groups"%(NumEdgeDiHydrox))
            NumSurfaceEpoxy = 0
        ## remove dihydrox, causing problems
        OxygenNumAndType["Edge_Hydroxyl"] += OxygenNumAndType["Dihydroxyl"]
        NumEdgeDiHydrox = 0
        
        ## Preferentially add hydroxyl groups to two-fold coordinated carbon atoms
        newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
        EdgeSites = recalculate_edge_sites(newcoords,abc,newatomtypes,False)
        write_pdb_file(newcoords,newatomtypes,abc,"edgecheck-2.pdb")
        print("Updated Edge Sites:")
        for Site in EdgeSites:
            print(Site,end=" ")
            try:
                PossibleDefectSites.remove(Site)
            except ValueError:
                pass
        print("")
        OxygenNumAndType["Edge_Hydroxyl"] += OxygenNumAndType["Surface_Hydroxyl"]
        NewO,NewH,EdgeSites,o_up,o_down = add_all_edgehydroxyl(NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
        #NewO,NewH,EdgeSites,o_up,o_down = add_all_hydroxyl_unrestricted(NewO,NewH,OxygenNumAndType,coords,abc,EdgeSites,NNList,o_up,o_down)
        NumEdgeHydrox = len(NewO) - NumEther - NumCarbAcid - NumCarbonyl - NumEdgeEpoxy - NumEdgeDiHydrox - NumSurfaceEpoxy
        print("Added %s Edge Hydrox groups"%(NumEdgeHydrox))
        newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
        write_pdb_file(newcoords,newatomtypes,abc,"edgecheck-3.pdb")

        ## Add remaining oxygen as surface hydroxyl groups
        OxygenNumAndType["Surface_Hydroxyl"] = OxygenNumAndType["Total"] - len(NewO)
        if OxygenNumAndType["Surface_Hydroxyl"] > 0:
            NewO,NewH,PossibleDefectSites,o_up,o_down = add_all_hydroxyl_unrestricted(NewO,NewH,OxygenNumAndType,coords,abc,PossibleDefectSites,NNList,o_up,o_down)
            NumSurfaceHydrox = len(NewO) - NumEther - NumCarbAcid - NumCarbonyl - NumEdgeEpoxy - NumEdgeDiHydrox - NumSurfaceEpoxy - NumEdgeHydrox
            print("Added %s Surface Hydrox groups"%(NumSurfaceHydrox))

    newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
    NewH = add_hydrogen(newcoords,newatomtypes,abc,NewH)
    newcoords,newatomtypes = add_new_atoms(coords,atomtypes,NewC,NewO,NewH)
    newcoords = wrap_coords(newcoords,abc)
    return(newcoords,newatomtypes)

######################################################################
## Add hydrogen to two-fold coordinated carbon
######################################################################
def add_only_hydrogen(coords,atomtypes,abc):
    NewH = add_hydrogen(coords,atomtypes,abc,[])
    coords,atomtypes = add_new_atoms(coords,atomtypes,[],[],NewH)
    newcoords = wrap_coords(coords,abc)
    return(coords,atomtypes)

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

def print_edges_for_vmd(coords,abc):
    NNList = nearest_neighbors_fv(coords,abc,2.9)
    EdgeSites = []
    for i in range(0,len(coords)):
        if len(NNList[i])< 12:
            EdgeSites.append(i)
    print("Edge site indices:")
    print("index ",end=" ")
    for i in EdgeSites:
        print(i,end=" ")
    print("")
    return()

######################################################################
## MAIN
######################################################################
def main():
    ## Read inputs 
    InputTileRelPath,OutputFile,PercentO,Fresh,OxygenProfile,BunchOrigin,BunchDiameter,DecorateEdgesFirst,ExcludeTileEdges = read_inputs()
    ## Reads input coords from pdb
    coords,atomtypes,abc = read_pdb_to_structure(InputTileRelPath)

    ## From % Oxygen and Hydroxyl:Epoxy ratio, determine amount of 
    ## hydroxyl and epoxy to place
    OriginalO,OriginalC,OriginalH = determine_oxygen(atomtypes)
    if Fresh:
        OxygenNumAndType = calculate_oxygen_num_type(OriginalC,OriginalO,PercentO,OxygenProfile)
    #OxygenNumAndType = manual_oxygen_num_type()
    elif Fresh == False:
        OxygenNumAndType = calculate_oxygen_num_type_aged(OriginalC,OriginalO,OxygenProfile,PercentO)
   
    ## Determine edge sites
    if DecorateEdgesFirst or Fresh != True:
        EdgeSites = define_defect_edge_sites(coords,atomtypes,abc)
        ## To select edge site indices in vmd
        #print_edges_for_vmd(coords,abc)

    ## Use all sites within a circle 
    if BunchDiameter != None:
        PossibleDefectSites = possible_sites_circle(coords,atomtypes,abc,BunchOrigin,BunchDiameter)
    ## use all sites
    else:
        if OriginalO == 0:
            PossibleDefectSites = use_all_defect_sites(atomtypes)
        else:
            PossibleDefectSites = recalculate_all_defect_sites(coords,atomtypes,abc)
        
    if ExcludeTileEdges:
        PossibleDefectSites = recalculate_defect_sites_excluding_edges(PossibleDefectSites,coords,atomtypes,abc)

    ## Place Oxygens and Hydrogens
    ## Fresh (hydroxyl and epoxy only) - for 1:1 epoxy:hydroxyl should work up to ~75%
    if Fresh:
        coords,atomtypes = add_oxygen_hydrogen(coords,atomtypes,abc,PossibleDefectSites,OxygenNumAndType)
    ## Aged (all func groups)
    elif Fresh == False and DecorateEdgesFirst != True:
        coords,atomtypes = add_oxygen_hydrogen_aged(coords,atomtypes,abc,PossibleDefectSites,EdgeSites,OxygenNumAndType)
    ## Aged (all func groups) -- use edge sites first
    elif Fresh == False and DecorateEdgesFirst == True:
        coords,atomtypes = add_oxygen_hydrogen_aged_edge(coords,atomtypes,abc,PossibleDefectSites,EdgeSites,OxygenNumAndType)
    ## Add hydrogen only
    else:
        coords,atomtypes = add_only_hydrogen(coords,atomtypes,abc)
        if PercentO != 0: 
            print("Warning: Percent oxygen to be added is nonzero, but neither --fresh nor --aged has been specified. Adding hydrogen only.")

    ## Write final tile geometry
    write_pdb_file(coords,atomtypes,abc,OutputFile) 

main()

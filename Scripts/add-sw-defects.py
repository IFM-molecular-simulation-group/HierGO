#!/usr/bin/python3
from argparse import ArgumentParser
import random
import math
import numpy as np

######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser(description='Add SW and GB defects to a graphene sheet.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--infile', nargs='?', help='input file name (pdb)',const='tile.pdb')
    #group.add_argument('--array', nargs='?', const=['20','20'], help='number of rows and \
    #columns of unit cells')
    parser.add_argument('--array', nargs=2, help='number of rows and \
    columns of tiles; if not set code defaults to 20,20')
    parser.add_argument('--sw', nargs='?', default=0.0, help='percent of SW defects',type=float) 
    parser.add_argument('--gb', nargs='?', default=0.0, help='percent of mini-GB defects',type=float)
    parser.add_argument('--exclude_defect_edges', action='store_true', help='optional \
    flag with no arguments to exclude hole edges as sites for topo defects')
    parser.add_argument('--outfile', nargs='?', default="topo-tile.pdb", type=str, help="output file name (pdb)")

    args = parser.parse_args()
    replicate = None
    if args.infile != None:
        print("Reading input structure from %s"%(args.infile))
    else:
        replicate = []
        if args.array == None:
            replicate.append(int(20))
            replicate.append(int(20))
        else:
            for i in args.array:
                replicate.append(int(i))
    return(args.infile,args.outfile,args.sw,args.gb,replicate,args.exclude_defect_edges)

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
def nearest_neighbor_point(i,coords,abc,cutoff):
    NNList = []
    ix = i[0]
    iy = i[1]
    iz = i[2]
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
## Create an orthographic box with a template sheet of C
######################################################################
def create_ortho_template(XRows,YRowsDouble):
    YRows = int(YRowsDouble/2)
    c1 = [ 1.841000,1.775000,0]
    c2 = [0.605000,2.484000,0]
    c3 = [0.605000,3.908000,0]
    c4 = [1.839000,4.620000,0]
    XRowCoords = []
    dx = 2.471
    for i in range(0,XRows):
        XRowCoords.append([float(c1[0]+dx*i),float(c1[1]),float(c1[2])])
        XRowCoords.append([float(c2[0]+dx*i),float(c2[1]),float(c2[2])])
        XRowCoords.append([float(c3[0]+dx*i),float(c3[1]),float(c3[2])])
        XRowCoords.append([float(c4[0]+dx*i),float(c4[1]),float(c4[2])])
    dy = 4.273000
    Coords = np.array([])
    AtomTypes = []
    for i in range(0,YRows):
        for b in range(0,len(XRowCoords)):
            if Coords.size == 0:
                Coords = np.append(Coords,[float(XRowCoords[b][0]),float(XRowCoords[b][1]+i*dy),float(XRowCoords[b][2])])
            else:
                Coords = np.vstack([Coords,[float(XRowCoords[b][0]),float(XRowCoords[b][1]+i*dy),float(XRowCoords[b][2])]])
            AtomTypes.append("C")
    #ABC = [float(YRowsDouble*2.13+1.234),float(XRows*2.46+2.471),63.350]
    ABC = [float(XRows*2.46),float(YRowsDouble*2.13),63.350]

    return(Coords,AtomTypes,ABC)
    
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
## Pick where to locate defect
######################################################################
def pick_defect_site(DefectTile):
    numatoms = len(DefectTile)
    DefectSite = random.choice(range(0,numatoms))
    #print("Defect Site:%s"%(DefectSite))
    neighbors = DefectTile.get_neighbors(site=DefectTile[DefectSite],r=1.45,include_image=True)
    # If the number of neighbors is less than 3 than the site chosen is 
    # adjacent to a defect, so another site is chosen
    if len(neighbors) < 3:
        DefectSite = random.choice(range(0,numatoms))
        neighbors = DefectTile.get_neighbors(site=DefectTile[DefectSite],r=1.45,include_image=True)
        if len(neighbors) < 3:
            sys.exit("Functionality to be added.")
    return(DefectSite)

######################################################################
## find the center of a series of sites
######################################################################
def td_center(Coords,Sites):
    origin = [ 0, 0, 0]
    for i in Sites:
        origin[0] += Coords[i][0]
        origin[1] += Coords[i][1]
        origin[2] += Coords[i][2]
    origin = [origin[0]/len(Sites),origin[1]/len(Sites),origin[2]/len(Sites)]
    return(origin)

######################################################################
## Check that the NN are within possible defect sites
######################################################################
def check_neighbors(DefectPossibleNums,Neighbors,Coords,ABC):
    XMax = ABC[0] - 1.8
    YMax = ABC[1] - 1.8
    PossibleNeighbors = []
    for i in Neighbors:
        if Coords[i][0] < XMax and Coords[i][1] < YMax and Coords[i][0] > 3 and Coords[i][1] > 3:
            PossibleNeighbors.append(i)
    return(PossibleNeighbors)

######################################################################
## 
######################################################################
def pick_sw_defect_site_ortho(Coords,AtomTypes,ABC,NNList,defect,DefectPossibleNums):
    tries = 0 
    check = False
    while check == False and tries < 100:
        # Choose initial defect site from possible sites
        DefectNum = random.choice(DefectPossibleNums)
        DefectNums = [DefectNum]
        # Choose neighboring defect site
        Neighbors = NNList[DefectNum]
        PossibleNeighbors = check_neighbors(DefectPossibleNums,Neighbors,Coords,ABC)
        if PossibleNeighbors != []:
            check = True
        tries += 1
    if check == False:
        print("Warning: not all defects placed.")
        return(DefectNums,DefectPossibleNumsNew,anchor)
    if defect == 2:
        DefectNum = random.choice(PossibleNeighbors)
        DefectNums.append(DefectNum)
        anchor = td_center(Coords,DefectNums)
    else:
        anchorassigned = False
        for i in Neighbors:
            if anchorassigned == False:
                dx = abs(Coords[i][0] - Coords[DefectNum][0])
                if dx < 1.0:
                    DefectNums.append(i)
                    anchor = td_center(Coords,DefectNums)
                    anchor[0] += 1.22
                    anchorassigned = True
        newneighbors = nearest_neighbor_point(anchor,Coords,ABC,3.80)
        DefectNums = []
        for i in newneighbors:
            DefectNums.append(i)
        DefectNums = DefectNums[0] #Why? Unknown

    DefectPossibleNumsNew = []
    for i in DefectPossibleNums:
        if i not in DefectNums:
            DefectPossibleNumsNew.append(i)
    return(DefectNums,DefectPossibleNumsNew,anchor)

######################################################################
## 
######################################################################
def create_sw_list(numAtoms,SWPercent):
    SWPossible = [ 24, 2 ]
    SWScaling = [ 41, 16 ]


    DefectList = []
    for i in range(0,len(SWPossible)):
        NumDefectAtoms = round(SWPercent[i]*numAtoms,0)
        NumDefects = int(round(NumDefectAtoms/SWScaling[i],0))
        for num in range(0,NumDefects):
            DefectList.append(SWPossible[i])
    return(DefectList)

######################################################################
## Rotate specified sites
######################################################################
def rotate_two_atoms(Coords,DefectNums,point):
    angle = np.pi/2

    print("Rotatating atoms:")
    for atom in DefectNums:
        print(atom,end=" ")
        px, py = Coords[atom][0], Coords[atom][1]
        ox, oy = point[0], point[1]
        qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
        qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
        Coords[atom] = [ qx, qy, point[2] ]
    print("")

    ## tweak coordinates for OPLS
    c1x = Coords[DefectNums[0]][0] 
    c1y = Coords[DefectNums[0]][1] 
    c2x = Coords[DefectNums[1]][0] 
    c2y = Coords[DefectNums[1]][1] 

    moved_x = c2x - c1x
    moved_y = c2y - c1y 
    angle = np.arctan(moved_y/moved_x)
    if moved_y == 0:
        bond_length = moved_x
    else:
        bond_length = moved_y / np.sin(angle)
    new_c2x = np.cos(angle) * (bond_length*0.8)
    new_c2y = np.sin(angle) * (bond_length*0.8)
    new_c1x = np.cos(angle) * (bond_length*0.2)
    new_c1y = np.sin(angle) * (bond_length*0.2)

    Coords[DefectNums[0]][0] = new_c1x + c1x
    Coords[DefectNums[0]][1] = new_c1y + c1y
    Coords[DefectNums[1]][0] = new_c2x + c1x
    Coords[DefectNums[1]][1] = new_c2y + c1y

    print("")
    print("about the anchor %s"%(point))
    return(Coords)

######################################################################
## Rotate specified sites
######################################################################
def rotate_atoms(Coords,DefectNums,point):
    theta = np.radians(90)
    # Rotation matrix
    r = np.array(( (np.cos(theta), -np.sin(theta)),
                   (np.sin(theta),  np.cos(theta)) ))
    print("Rotatating atoms:")
    for atom in DefectNums:
        print(atom,end=" ")
        print(Coords[atom])
        v = np.array((Coords[atom][0]-point[0], Coords[atom][1]-point[1]))
        v_r = r.dot(v)
        Coords[atom] = [ v_r[0]+point[0], v_r[1]+point[1], point[2] ]
    print("")
    print("about the anchor %s"%(point))
    return(Coords)

######################################################################
## Only include sites within a circle of radius 50 Ang within center
######################################################################
def possible_sites_ortho_circle(Tile,XRows,YRows):
    anchorx = XRows*2.46/2
    anchory = YRows*2.13/2
    anchor = [anchorx,anchory,0]
    #anchor = [51.879002, 46.756001, 31.675]
    PossibleSites = Tile.get_sites_in_sphere(pt=anchor,r=25,include_image=False,include_index=True)
    DefectPossibleNums = []
    DefectPossibleNumsLarge = []
    for i in PossibleSites:
        ## To remove sites border holes
        #neighbors = Tile.get_neighbors(site=Tile[i[2]],r=2.86,include_image=True,include_index=True)
        #if len(neighbors) >= 12:
            #DefectPossibleNums.append(i[2])
        #lessen liklihood of 24 atom defect being placed strangely
        newneighbors = Tile.get_neighbors(site=Tile[i[2]],r=3.80,include_image=True,include_index=True)
        if len(newneighbors) >= 18:
            DefectPossibleNumsLarge.append(i[2])
        else:
            neighbors = Tile.get_neighbors(site=Tile[i[2]],r=2.86,include_image=True,include_index=True)
            if len(neighbors) >= 12:
                DefectPossibleNums.append(i[2])
        ## all sites
        #DefectPossibleNums.append(i[2])
    print("Possible sites for 2 atom rotation:")
    for atom in DefectPossibleNums:
        print(atom,end=" ")
    print("")
    print(DefectPossibleNums)
    return(DefectPossibleNums,DefectPossibleNumsLarge)

######################################################################
## All sites are possible sites except outermost
######################################################################
def possible_sites_no_outer(Coords,ABC,cutoff,ExcludeDefectEdges):
    if ExcludeDefectEdges:
        NNList = nearest_neighbors_fv(Coords,ABC,1.44)

    XMax = ABC[0] - cutoff
    YMax = ABC[1] - cutoff
    DefectPossibleNums = []
    HoleEdges = []
    count = 0
    for atom in range(0,len(Coords)):
        i = Coords[atom]
        if i[0] < XMax and i[1] < YMax and i[0] > 3 and i[1] > 3:
            if ExcludeDefectEdges:
                if len(NNList[atom]) > 2:
                    DefectPossibleNums.append(count)
                else:
                    HoleEdges.append(count)
                    for NN in NNList[atom]:
                        HoleEdges.append(NN)
                        ## exclude NNN
                        for NNN in NNList[NN]:
                            HoleEdges.append(NNN)
            else:
                DefectPossibleNums.append(count)
        count += 1
    if ExcludeDefectEdges:
        NewDefectPossibleNums = []
        for i in DefectPossibleNums:
            if i not in HoleEdges:
                NewDefectPossibleNums.append(i)
        return(NewDefectPossibleNums)

    ## Print possible sites
    #print("Possible defect sites:")
    #for i in DefectPossibleNums:
    #    print(i,end=" ")
    #print("")
    return(DefectPossibleNums)

######################################################################
## TODO: Make rings nice for OPLS
######################################################################
def fix_up_rings(Coords,ABC,NNList):
    NNList = nearest_neighbors_fv(Coords,ABC,1.44) ##new NN list 
    for i in range(0,len(Coords)):
        if len(NNList[i]) > 4: 
            print(i)
            #do something
    return(Coords)

######################################################################
## Write an output PDB file
######################################################################
def write_pdb_file(Coords,AtomTypes,ABC,FileName):
    newlines = []
    count = 0
    for i in range(0,len(Coords)):
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,AtomTypes[i],"UNL",1,Coords[i][0],Coords[i][1],Coords[i][2],1.00,0.00,AtomTypes[i]))
        count += 1

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
    InputFile,OutputFile,SWPercent,GBPercent,Array,ExcludeDefectEdges =  read_inputs()
    if InputFile == None:
        Coords,AtomTypes,ABC = create_ortho_template(Array[0],Array[1])
        write_pdb_file(Coords,AtomTypes,ABC,"graphene.pdb")
    else:
        Coords,AtomTypes,ABC = read_pdb_to_structure(InputFile)

    SWDefectList = create_sw_list(len(AtomTypes),[GBPercent/100.0,SWPercent/100.0])
    DefectPossibleNums = possible_sites_no_outer(Coords,ABC,3.2,ExcludeDefectEdges)

    count = 0
    string = ""
    ## DefectNums are the atom indexs to rotate
    ## DefectPossibleNums are the atoms which can be selected as defect sites
    for defect in SWDefectList:
        #DefectNums,origin = pick_sw_defect_site(Tile,defect,DefectPossibleNums)
        NNList = nearest_neighbors_fv(Coords,ABC,1.44)
        if defect == 2:
            DefectNums,DefectPossibleNums,origin = pick_sw_defect_site_ortho(Coords,AtomTypes,ABC,NNList,defect,DefectPossibleNums)
            Coords = rotate_two_atoms(Coords,DefectNums,origin)
        elif defect == 24:
            DefectPossibleNums = possible_sites_no_outer(Coords,ABC,8.0,ExcludeDefectEdges)
            DefectNums,DefectPossibleNums,origin = pick_sw_defect_site_ortho(Coords,AtomTypes,ABC,NNList,defect,DefectPossibleNums)
            Coords = rotate_atoms(Coords,DefectNums,origin)
        write_pdb_file(Coords,AtomTypes,ABC,"defect-%s.pdb"%(count)) #de-bugging
        count +=1

    ## TODO
    #Coords = fix_up_rings(Coords,ABC)

    #print(string)
    write_pdb_file(Coords,AtomTypes,ABC,OutputFile) 
    #write_output_file(Coords,AtomTypes,ABC,"defect-final.cif")

main()

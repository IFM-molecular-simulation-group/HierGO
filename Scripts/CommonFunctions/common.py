#!/usr/bin/python

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
## FV NN
######################################################################
def distance(x0, x1, box, pbc):
    # xo is a position of one atom, x1 is an array of positions
    # use the pbc bool mask to set the periodicity
    delta = np.abs(x0 - x1)
    delta[:,pbc] -= box[pbc] * np.round(delta[:,pbc]/(box[pbc]))
    return np.sqrt((delta ** 2).sum(axis=-1))

def nearest_neighbors_fv(coords, box, cutoff):
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
## Generate nearest neighbor list from coords
######################################################################
def nearest_neighbor_single(i,coords,abc,cutoff):
    ix,iy,iz = i
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
## Print edges for vmd visualizations
######################################################################
def print_edges_for_vmd(coords,abc):
    NNList = nearest_neighbors(coords,abc,2.9)
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

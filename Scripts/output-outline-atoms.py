#!/usr/bin/python3
from argparse import ArgumentParser
## Find the defect outline atoms, ignore oxygen


######################################################################
## Read inputs
######################################################################
def read_inputs():
    parser = ArgumentParser()
    parser.add_argument('-infile', type=str)

    args = parser.parse_args()
    return(args.infile)

######################################################################
## Read PDB
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
## Define defect edge sites
######################################################################
def define_defect_edge_sites(coords,atomtypes,abc):
    DefectEdgeSites = []
    for i in range(0,len(coords)):
        if atomtypes[i] == "C1":
            ix = coords[i][0]
            iy = coords[i][1]
            iz = coords[i][2]
            neighbors = 0
            for j in range(0,len(coords)):
                if atomtypes[j] == "C1":
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
                    if d <= 2.86:
                        neighbors += 1
            if neighbors < 13:
                #print("C          {:8.3f}{:8.3f}{:8.3f}".format(Tile[i].coords[0],Tile[i].coords[1],0.00)) #de-bugging
                DefectEdgeSites.append(i)
    return(DefectEdgeSites)

######################################################################
## Write an output PDB file
######################################################################
def write_pdb_file(InputStructure,FileName):
    newlines = []
    count = 0
    for i in InputStructure:
        if "C" in str(i.species):
            i.coords[2] = 0 #z always in plane for Cs
        newlines.append("{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM",count,str(i.species),"UNL",1,i.coords[0],i.coords[1],i.coords[2],1.00,0.00,str(i.species)))
        count += 1

    file = open(FileName,"w")
    file.write("AUTHOR    GENERATED BY NG TILE DECORATION\n")
    file.write("CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P1          1\n".format(InputStructure.lattice.abc[1],InputStructure.lattice.abc[2],InputStructure.lattice.abc[0],InputStructure.lattice.beta,InputStructure.lattice.gamma,InputStructure.lattice.alpha))
    for line in newlines:
        file.write("%s\n"%(line))
    file.close()

    return()

######################################################################
## MAIN
######################################################################
def main():
    filename = read_inputs()
    coords,atomtypes,abc = read_pdb_to_structure(filename)
    DefectEdgeSites = define_defect_edge_sites(coords,atomtypes,abc)
    print("")
    for i in DefectEdgeSites:
        print(i,end=" ")
    #write_pdb_file(struct,"output.pdb")
main()

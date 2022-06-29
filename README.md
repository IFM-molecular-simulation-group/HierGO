# HierGO 

## Getting started
```
## This workflow requires git. 
## Getting started with GIT: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

## Inputs to all scripts can be reviewed using the -h option:
python3 ../Scripts/create-tile.py --h
```

## Creating a GO tile with equal epoxy:hydroxyl ratio
```
## Create a tile 
python3 ../Scripts/create-tile.py --percentvac 10
## Decorate tile
python3 ../Scripts/decorate-tile.py --fresh --percento 30
```

## Creating a GO tile with all functional groups and topo defects
```
## Create a tile
python3 ../Scripts/create-tile.py --percentvac 20
## Add topo defects 
python3 ../Scripts/add-sw-defects.py --sw 5 --gb 5
## Decorate tile
python3 ../Scripts/decorate-tile.py --infile topo-tile.pdb --aged --percento 20
```

## Stitch tiles together
```
python3 ../Scripts/stitch.py --infiles file1 file2 file3 file4
## Manually check the geometry for questionable atoms, if they occur remove them
python3 ../Scripts/remove-atoms.py -i 1 2 3 4 5 
```

## Atom typing for OPLS
```
## Run the atom typing script for use with OPLS
##   **Unfortunately, this script is quite slow for a large number of atoms, feel free to look into 
##     making it faster
python3 ../Scripts/atom-typing.py --infile removed.pdb
## If there are unbound atoms flagged by the typing script, 
## check where they are, and if they should not be present remove them
## Ideally, you would have removed or adjusted these after your manual 
## check of the stitched geometry
python3 ../Scripts/remove-atoms.py -i 1 2 3 4 5
## Re-run typing script until no unbound atoms occur
python3 ../Scripts/atom-typing.py --infile removed.pdb
```

## Directory layout
```
Directory          Description
---------------------------------------------------------------------------------------------
HoleLibrary/       Current library of .xyz holes
Scripts/           All scripts
```

# HierGO 

## Quick guide
---
```
## this workflow requires git. 
## getting started with git: https://git-scm.com/book/en/v2/getting-started-installing-git

## Inputs to all scripts can be reviewed using the -h option:
python3 ../Scripts/create-tile.py --h
```
### Creating a GO tile with equal epoxy:hydroxyl ratio
```
## Create a tile 
python3 ../Scripts/create-tile.py --percentvac 10
## Decorate tile
python3 ../Scripts/decorate-tile.py --fresh --percento 30
```

### Creating a GO tile with all functional groups and topo defects
```
## Create a tile
python3 ../Scripts/create-tile.py --percentvac 20
## Add topo defects 
python3 ../Scripts/add-sw-defects.py --sw 5 --gb 5
## Decorate tile
python3 ../Scripts/decorate-tile.py --infile topo-tile.pdb --aged --percento 20
```

### Stitching tiles together
```
python3 ../Scripts/stitch.py --infiles file1 file2 file3 file4
## Manually check the geometry for questionable atoms, if they occur remove them
python3 ../Scripts/remove-atoms.py -i 1 2 3 4 5 
```

### Atom typing for OPLS
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

## User manual
---
### Step 1: Create a tile
```
python3 ../Scripts/create-tile.py
```
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --percentvac | A positive number or zero | The percentage of the vacancy (default is 5 percent). | --percentvac 10 |
| --distribution | small, medium, large or mixed | The size of the hole on the tile (default is mixed). | --distribution mixed |
| --outfile | Filename | Naming the output of the file (default is tile.pdb) | --outfile tile\_10.pdb |
| --xfactor | A positive number | Number of replicate of the unit cell in x direction (default is 20) | --xfactor 25 |
| --yfactor | A positive number | Number of replicate of the unit cell in y direction (default is 20) | --yfactor 25 |
| --keep_hanging_carbon | None | Keep hanging carbon items, if not specified hanging carbon items will be trimmed | N/A |
| --non_periodic | None | Each tile only can be used individually and cannot be stitches together directly. | N/A |

*Warning:* 
The actual amount of vacancy of the tile may not be the same as the specified amount of vacancy due to the trimming of the carbons.

### Step 2: Add defect
`python3 ../Scripts/add-sw-defects.py`

| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --infile | Filename from last step | Filename from last step | --infile tile\_10.pdb |
| --array | Two positive integers | Create a tile with zero vacancy | --array 20 20 |
| --sw | A positive number or zero | Percent of SW defects (default is 0) | --sw 5 |
| --gb | A positive number or zero | Percent of mini-GB defects (default is 0) | --gb 4 |
| --exclude\_defect\_edges| None | Exclude vacancy edges as sites for topological defects | N/A |
| --outfile | Filename | Naming the output of the file (default is topo-tile.pdb) | --outfile tile\_topo.pdb |

*Warning:* 
If multiple tiles need to be stitched together, this code needs to run even if the percentage of defects in a given tile is zero. This ensures that the tiles are created with the same dimensions.

## Step 3: Decorate oxygen
`python3 ../Scripts/decorate-tile.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --infile | Filename from last step | Filename from last step | --infile tile\_topo.pdb |
| --percento | Positive number or zero | Percentage of oxygen (default is 0) | --percento 10 |
| --fresh | A ratio | Decorate with hydroxyl and epoxy groups only in a ratio specified by hydroxyl:epoxy (default is 1:1) | --fresh 2:1 |
| --aged | 600, 900, 1200 or 2500 | Decorate with hydroxyl, epoxy, carboxylic acid, epoxide and non-epoxide ether in a ratio specified by annealing temperature (default is 900). | --aged 900 |
| --outfile | Filename | Naming the output of the file (default is dec-tile.pdb) | --outfile tile\_dec.pdb |
| --decorate\_edges\_first | None | Decorate defect edges before the basal plane | N/A |
| --bunching | A positive number or zero plus a positive number | The first number indicates the index of the atom. The second number indicates the diameter of circle in Angstroms that where the oxygen is going to be decorated. | --bunching 309 10 |
| --exclude\_tile\_edges| None | Exclude tile edges in the decoration process | N/A |

*Warnings:* 
1. If multiple tiles need to be stitched together, this code needs to run even if the percentage of oxygen in a given tile is zero. This ensures that the tiles are created with the same dimensions.
2. Fresh or aged must be specified in the options. Otherwise, no oxygen will be added.
3. The code may fail when adding higher amount of oxygen (such as 40%). Please try a few times.
4. For the bunching option, the code cannot produce the correct result when the diameter is too small. Please choose a larger diameter.

## Step 4: Stitch
`python3 ../Scripts/stich.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --infiles | Filenames | Filenames that stitch together | --infile tile\_10.pdb tile\_2.pdb tile\_5.pdb |
| --percento | Positive number or zero | Oxygen content of the product (default is 0) | --percento 15 |
| --array | Two positive numbers | The dimension of array [x,y]. If this flag is not provided, the code will attempt to make a square array given a proper number of tiles have been inputted. | --array 3 1 |
| --keep\_overlap | None | Keep overlapping oxygen and carbon (by default these atoms are trimmed) | N/A |
| --add\_oxygen\_to\_sheet\_edges | None | Add oxygen determined from the --percento flag to the sheet edges only (by default is false) | N/A |
| --keep\_hanging\_carbon | None | Keep hanging carbon atoms (by default these atoms are trimmed) | N/A |
| --ordered\_tiles | None | Keep the input order of the tiles | N/A |
| --pbc | None | Keep periodic boundary conditions, otherwise buffer distance of 10 Angstroms is added | N/A |
| --outfile | Filename | Naming the output of the file (default is stich.pdb) | --outfile tile\_stitch.pdb |

*Warnings:*
1. When adding oxygen to the product, the oxygen will be located at the edges first.
2. The code may fail if you trying to make too great adjustment to the total oxygen contents. For example, an adjustment of 1-2% should be fine.
3. Manually check the geometry for questionable atoms. If they occur, remove them.

## Step 5: Remove atoms
`python3 ../Scripts/remove-atoms.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --infile | Filename | Filename that have atom(s) that need to be removed | --infile tile\_10.pdb |
| --i | Positive numbers or zero | The index of atom(s) that need to be removed | --i 7 8 90 |
| --outfile | Filename | Naming the output of the file (default is removed.pdb) | --outfile tile\_remove.pdb |

*Warning:*
This step can be done after conducting any step above. It does not have to be at the last step.

## Step 6: Atom typing
`python3 ../Scripts/atom-typing.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --infile | Filename | Name of the input file | --infile tile\_10.pdb |
| --outfile | Filename | Name the output of the file (default is typed.pdb) | --outfile tile\_atom.pdb |
| --pbc | None | Default to ignoring periodic boundary condition. Only required when a single tile is made to be used. Buffer distance is 10 Angstroms for pbc. | N/A |

*Warning:*
Some warning could appear especially for the stitched sheet. It is better to check possible problem atoms in the visualization software.


## Optional routines:
## Rupture
`python3 ../Scripts/create-rupture.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --size | A positive number | Number of replicate of the unit cell in both x and y direction (default is 20) | --size 40 |
| --outfile | Filename | Name the output of the file (default is rupture.pdb) | --outfile tile\_rup.pdb |
| --atoms | Positive number or zero | Number of oxygen in the product (default is 4). | --atoms 6 |

*Warnings:*
1. The generated tile has zero vacancy.
2. If you want to decorate oxygen using the rupture tile from this code, it could produce wrong structure when adding a higher amount of oxygen (such as 35%).

## Trim hanging carbons
`python3 ../Scripts/trim-hanging-Cs.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| --infile | Filename | Name of the input file | --infile tile.pdb |
| --outfile | Filename | Name of the output file | --outfile tile\_rup.pdb |
| --atoms | Positive number or zero | Number of atoms need to be removed | --atoms 6 |
| --orthogonalize | None | Change the periodic boundary to orthogonal | N/A |
*Warnings:*
1. The number of atoms need to be provided otherwise the code cannot work appropriately.
2. There is a default number of atoms to be removed from a structure. If the number of the atoms is below that number, the default number will be applied. Otherwise the provided number will be applied (the actual number of atoms removed could be slightly different to the provided number due to the intrinsic structure)

## Translate pdb file according by a given vector
`python3 ../Scripts/move\_pdb\_by\_vector.py`
| **Option names** | **Input value** | **Meaning** | **Example** |
| --- | --- | --- | --- |
| -h or -help | None | Show the help message | N/A |
| -i | Filename | Name of the input file | -i repture.pdb |
| -out | Filename | Name of the output file | -out tile\_mv.pdb |
| -vector | A vector | Move the pdb file based on the vector | -vector [10,7,9] |

*Warning:*
The cell must be orthogonal.

## Directory layout
```
Directory          Description
---------------------------------------------------------------------------------------------
HoleLibrary/       Current library of .xyz holes
Scripts/           All scripts
```

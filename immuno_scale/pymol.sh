## bin/batch/pymol

### fetch biological unit of pre-fustion F

fetch 5UDE, type=pdb1
split_state 5UDE

## display surface
show surface
# color the surface white
set surface_color, gray70, *


### color epitope residue


color firebrick, resi 33-41 
color firebrick, resi 44-52
color limegreen, resi 53-61, 
color limegreen, resi 57-65, 
color limegreen, resi 86-94
color skyblue, resi 140-152
color yelloworange, resi 164-172, 
color yelloworange, resi 178-187, 
color yelloworange, resi 203-211, 
color yelloworange, resi 221-231, 
color yelloworange, resi 233-241, 
color yelloworange, resi 250-260, 
color yelloworange, resi 271-288, 
color yelloworange, resi 299-307
color violet, resi 314-322
color teal, resi 407-415, 
color teal, resi 431-439
color orange, resi 499-507


## postfusion
fetch 3RRR, type=pdb1
split_state 3RRR
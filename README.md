# ambitus-postprocessor
This is a new version of the postprocessing script included in the original ambitus for geographic range reconstruction. The previous version was hard-coded for two-node comparisons; this new version loops through all climate scenarios and all nodes in the phylogeny. It will not take into account the age of nodes in the phylogeny, so keep what is relevant and remove the rest.

Usage example, where the ambitus run is in the parent directory: 
```
./summarize_nodal_distributions.py 6 ./../ambitus_run 24 --list cclgmbi ccmidbi lig_2_5min_bio_
```

The first argument is the number of niche variables (numbered 1-6 following the naming conventions of `ambitus`. The second is a path to a completed `ambitus` run. The third is the number of nodes in the phylogeny (which can be fetched from `ambitus` output or by counting entries in the file `reconstructionnodelist.command` that is produced when the program is run. The last is a list of climate scenario files in the `./scenarios/` subdirectory of the working directory. The file names should be in three parts like so: 'arbitrary prefix' 'variable#' '.tif'. The list pertains to the prefixes.

The program will then perform all the calculations in succession. In a tree with 24 nodes, 6 variables, and 3 scenarios, this program took about 36 hours. If there is interest, I could rewrite this to accept a list of nodes rather than reconstructing on all of them.

Prerequisites are R package raster and GDAL (specifically the executables).

The main prototype of the code is done but there are some minor issues. 


# TODO 
1) Distribute wfcs across all the kpools (now data is duplicated)
2) Rearrange non-local part (compute beta functions at the start and blocking (loop tiling))
3) Prepare test cases
4) Improve OPENMP


# Some Uncertain improvements 
1) dVloc routines needs to be slightly rewritten to avoid multiple array allocations


# Known issues
1) Currently, the code assumes presence of SOC for non-collinear calculations.
    which implies that f-coeffs must be computed for soc case only 



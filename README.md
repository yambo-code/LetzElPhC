The main prototype of the code is done but there are some minor issues. 
The developement will be very slow from now on. I will only develop in 
my free times.

# TODO (Short term). (Should be addressed after return)
1) Distribute wfcs across all the kpools (now data is duplicated)
2) Compute the representation matrices to rotate in q space
3) Prepare test cases

# TODOs (LONG TERM)
1) Prepare an interface 
2) read dvscf directlty from binary


# Some Uncertain improvements (should be addressed after return)
1) Improve IO performace(AION and IRIS FS report very bad file open times)
2) fft transpose.c routine ( FFTW MPI has a direct inplace transpose but looks 
like there is a bug if the Block size is not set to default). Alternatively,
if possible, directly perform transpose without rearrangement ?
3) dVloc routines needs to be slightly rewritten to avoid multiple array allocations




################################################################################
#                                                                              #
#                                                                              #
#            HORSES3D High-Order (DG) Spectral Element Solver                  #
#                                                                              #
#                                                                              #
################################################################################

To build and run:

1. cd to any of the test cases in this directory
2. compile the ProblemFile.f90 in SETUP
	make "Compilation options used to compile main code" (NOTE MKL is required for some of the test cases)
3. execute
	./horses3d.mu NAMEOFCONTROLFILE.control
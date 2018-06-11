**** TREx build and run information ****

- Please ensure correct versions of modules are being used, if you are using a Warwick Machine these can be accessed via a few lines of code that can be sourced from a .sh file: 

	module use --append /software/epp/modules
	module load brews/epp
 
	export LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}


// Building

- TREx is built using the 'make' command simply type this into the terminal in the src directory.


// Running

- First input data needs to be voxelized. This is done using the Voxelize_CCD_XY.cc script.

- Load root and input the following commands:

	.L macros/Voxelize_CCD_XY.cc+
	Voxelize(<"data file">,xdim,ydim,0,<"Ideal" or "Real">) 

  Whereby:
  xdim = voxel dimention of x axis in units e^-5 m
  ydim = voxel dimention of y axis in units e^-5 m

- Next TREx is run using: 

	./RunTRExReconCCD.exe <voxelized file> 


// Running on the cluster 

- Sammy has created a simple bash script to help with running TREx on the Warwick cluster, this can be found in runTREx.sh and can be run using:

	source runTREx.sh <input files>
  
  this script includes both voxelisation and running of TREx.






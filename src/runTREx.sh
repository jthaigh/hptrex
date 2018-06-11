#!/bin/sh

# INFO: Default setting for voxelize is ("<inputfile>",234,234,0,"<Ideal or Real>")
# INFO: Designed to run in the hptrex/src working directory, paths will need to be editted to run elsewhere
# INFO: How to run: source runTREx.sh <location of data> 

# Clear any previous voxelized files and create a temporary location for more
rm -rf voxelized/
mkdir voxelized

# Voxelizing

for file in $*
do

echo "Voxelising file: "$file
bnfile=$(basename "$file")
echo "Basename of file: "$bnfile
exbnfile=${file%.root} 
echo "Extracted basename of file: "$exbnfile


root -l << EOF
.L macros/Voxelize_CCD_XY.cc++
Voxelize("$file",277,277,0,"Ideal")     
.q
EOF

mv $exbnfile"_voxelsIdeal_277e-5m_277e-5m.root" voxelized/
done


echo -------- Voxelized --------



# Running TREx in submitted jobs

voxeldir=voxelized/*

for f in $voxeldir
do
echo "Processing file: "$f 
bsub -q xlong -e "error`basename $f | sed 's/\..*//'`.txt"  -J "`basename $f | sed 's/\..*//'`" ./RunTRExReconCCD.exe $f
done



echo -------- Finished --------




#!/bin/sh

TRExdir=$PWD

filepath=$1

outputpath=$2

fbname=$(basename "$1")

echo copying $fbname

cp $filepath $outputpath

name="Recon_$fbname"

echo Renaming File to: $name

mv $outputpath/$fbname $outputpath/$name  

echo Running TREx Reconstruction

$TRExdir/RunTRExRecon.exe $outputpath/$name
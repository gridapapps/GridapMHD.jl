#!/bin/bash

#Use Julia to compute the ratios for the BL geomtric clusterings

julia -e\
'
using DelimitedFiles

include("pass_params.jl")	#The parameters are read from this file
include("../../meshes/tube_block/params_calculator.jl")

ratio=mesh_params(Ha = Ha, N = N_r, n=n, a=(R-p))

pass=(string("r=",ratio)*";",
	  string("R=",R)*";",
	  string("N_r=",N_r)*";",
	  string("n=",n)*";",
	  string("N_a=",N_a)*";",
	  string("N_L=",N_a)*";",
	  string("L=",L)*";",
	  string("p=",p)*";"
	  )

writedlm("mesh_params",pass)

'

#Edit the gmsh input file with the calculated information

#Eliminate old information (if any)
sed -i '/r=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/R=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/N_r=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/n=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/N_a=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/N_L=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/L=/d' $mesh_path/tubeMaPLE-block.geo
sed -i '/p=/d' $mesh_path/tubeMaPLE-block.geo

#Write the new information

sed -e '/Mesh parameters/r./mesh_params' $mesh_path/tubeMaPLE-block.geo > $mesh_path/aux.geo
mv $mesh_path/aux.geo $mesh_path/tubeMaPLE-block.geo

rm mesh_params 

#Generate the mesh using gmsh

gmsh $mesh_path/tubeMaPLE-block.geo -3 -o $mesh_path/../tube_MaPLE.msh

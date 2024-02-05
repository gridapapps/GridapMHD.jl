#!/bin/bash

export path=$GRIDAPMHD/meshes/expansion

#Use Julia to compute the ratios for the BL geomtric clusterings

julia -e\
'
using DelimitedFiles

include("input_params.jl")
include("../../meshes/expansion/params_calculator.jl")
params=mesh_params(Ha = Ha, N_Ha = N_Ha, n_Ha=n_Ha, N_side=N_side, n_side=n_side,β=1/r)

pass=(string("R_Ha=",params[1])*";",
	  string("R_side=",params[2])*";",
	  string("N_Ha=",N_Ha)*";",
	  string("n_Ha=",n_Ha)*";",
	  string("N_side=",N_side)*";",
	  string("n_side=",n_side)*";",
	  string("n_inlet=",n_inlet)*";",
	  string("n_outlet=",n_outlet)*";",
	  string("r_exp=",R_exp)*";",
	  string("r=",r)*";",
	  string("betta=",betta)*";",
	  string("L_in=",L_in)*";",
	  string("L_out=",L_out)*";"
	  )

writedlm("mesh_params",pass)
'

#Edit the gmsh input file with the calculated information

#Eliminate old information
sed -i '/R_Ha=/d' $path/Expansion_B.geo
sed -i '/R_side=/d' $path/Expansion_B.geo
sed -i '/N_Ha=/d' $path/Expansion_B.geo
sed -i '/N_side=/d' $path/Expansion_B.geo
sed -i '/n_Ha=/d' $path/Expansion_B.geo
sed -i '/n_side=/d' $path/Expansion_B.geo
sed -i '/n_inlet=/d' $path/Expansion_B.geo
sed -i '/n_outlet=/d' $path/Expansion_B.geo
sed -i '/r_exp=/d' $path/Expansion_B.geo
sed -i '/r=/d' $path/Expansion_B.geo
sed -i '/betta=/d' $path/Expansion_B.geo
sed -i '/L_in=/d' $path/Expansion_B.geo
sed -i '/L_out=/d' $path/Expansion_B.geo

#Write the new information

sed -e '/Mesh parameters/r./mesh_params' $path/Expansion_B.geo > $path/Expansion_B-2.geo
rm mesh_params
mv $path/Expansion_B-2.geo $path/Expansion_B.geo

#Generate the mesh

gmsh $path/Expansion_B.geo -3 -o $path/Expansion_computed.msh



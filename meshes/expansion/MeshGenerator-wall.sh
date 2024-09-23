#!/bin/bash

export path=$GRIDAPMHD/meshes/expansion

julia -e\
'
using DelimitedFiles

include("input_params.jl")
include("params_calculator.jl")
R_1= R_calc(;N=N_1,Dx=dx,L=L_out)
R_2= R_calc(;N=N_2,Dx=dx,L=L_in)

pass=(string("Ha=",Ha)*";",
        string("N_Ha=",N_Ha)*";",
        string("N_s=",N_s)*";",
        string("exp_Ha=",exp_Ha)*";",
        string("exp_s=",exp_s)*";",
        string("r=",r)*";",
        string("betta=",betta)*";",
        string("L_in=",L_in)*";",
        string("L_out=",L_out)*";",
        string("t=",t)*";",
        string("N_1=",N_1)*";",
        string("N_2=",N_2)*";",
        string("N_w=",N_w)*";",
        string("R_1=",R_1)*";",
        string("R_2=",R_2)*";"
        )

writedlm("mesh_params",pass)
'

#Edit the gmsh input file with the calculated information



#Eliminate old information
sed -i '/Ha=/d' Expansion_B-wall.geo
sed -i '/N_Ha=/d' Expansion_B-wall.geo
sed -i '/N_s=/d' Expansion_B-wall.geo
sed -i '/exp_Ha=/d' Expansion_B-wall.geo
sed -i '/exp_s=/d' Expansion_B-wall.geo
sed -i '/N_1=/d' Expansion_B-wall.geo
sed -i '/N_2=/d' Expansion_B-wall.geo
sed -i '/N_w=/d' Expansion_B-wall.geo
sed -i '/R_1=/d' Expansion_B-wall.geo
sed -i '/R_2=/d' Expansion_B-wall.geo
sed -i '/r=/d'   Expansion_B-wall.geo
sed -i '/betta=/d' Expansion_B-wall.geo
sed -i '/L_in=/d' Expansion_B-wall.geo
sed -i '/L_out=/d' Expansion_B-wall.geo
sed -i '/t=/d' Expansion_B-wall.geo

#Write the new information

sed -e '/Mesh parameters/r./mesh_params' Expansion_B-wall.geo > Expansion_B-2.geo
mv Expansion_B-2.geo  Expansion_B-wall.geo

rm mesh_params
rm input_params.jl

#Generate the mesh

gmsh $path/Expansion_B-wall.geo -3 -o $path/../Expansion-wall_computed.msh



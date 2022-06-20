# Abrupt expansion along the direction of the magnetic field
# Validation problem with respect to the experiment published 

using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapMHD:main

#----------Problem setting-------


function ____expansion(;
   Ha = 1.0,                            #Hartmann number
   Re = 1.0,                            #Reynolds number
   c_w = 0.0727,                        #Wall conductance ratio
   B_dir = VectorValue(0.0,1.0,0.0),    #Direction of the B-field (norm =1)
   α_w = 100,                           #Penalty for the wall boundary condition
   a = 1,                               #Characteristic Length
   L = 6,                               #Half length of the channel (normalized by a)
   ρ = 1,                               #Density (Kg/m3)
   ν = 1,                               #Kinematic viscosity (Pa·s)
   σ = 1,                               #Electrical conductivity (S/m)
   np = 2                               #Number of parts for parallel computations
)

# Derived variables

  N = Ha^2/Re                      #Interaction parameter or Stuart number

  U0 = Re*ν/R                      #Velocity scale
  B0 = Ha*sqrt((ρ*ν)/σ)/R          #Module of the external B-field


#------------Mesh and geometry------------------------

  #This is for parallel id
 ## parts = get_part_ids(sequential,1)
  parts = np 

  t = PTimer(parts,verbose=true)

  model = GmshDiscreteModel("../geometries/Expansion_coarse.msh") 
  Ω = Triangulation(model)
  

#-------Dictionary of parameters for the solver-----------

  # Parabolic profile at the inlet
  function g_u(coord,a,L)
     r = sqrt(coord[1]^2 + coord[2]^2)
     zz = coord[3]
     tol = 1e-6
     if zz < (-L/a + tol)
         return VectorValue(0.0, 0.0, U0)
     else
         return VectorValue(0.0, 0.0, 0.0)
     end
  end

  #Thin wall boundaries

  Γ_w   = Boundary(model, tags="wall")

  z = zero(VectorValue{3,Float64})

       params = Dict(
          :ptimer=>t,
	  :debug=>true,
          :fluid=>Dict(
            :domain=>model,
            :α=>1,
            :β=>1,
            :γ=>1,
            :u=>Dict(
              :tags=>["wall","inlet"],
              :values=>[z, g_u]),
            :j=>Dict(
              :tags=>["inlet","outlet","wall"],
              :values=>[z, z, z]),
            :f=>z,
            :B=>VectorValue(0.0,1.0,0.0),
          ),
        )

#-----------------Solver call-------------------------

xh = main(params)

#----------------Writing results--------------------

uh, ph, jh, φh = xh

writevtk(Ω,"Expansion_debug",cellfields=["u"=>uh,"p"=>ph,"j"=>jh,"phi"=>φh])


